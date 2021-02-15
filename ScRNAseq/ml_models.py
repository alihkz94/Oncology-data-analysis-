#!/usr/bin/env python
"""
Train a machine learning model to predict treatment response based on scRNA-seq transcriptional signatures.
This script uses scVI to obtain latent representations from single-cell data and trains a PyTorch
MLP classifier to perform binary classification (e.g., responder vs. non-responder).

Pipeline Overview:
1. Load a processed AnnData (.h5ad) file containing scRNA-seq data and a binary 'treatment_response' label in adata.obs.
2. Set up and train a scVI model to learn a low-dimensional latent space.
3. Extract the latent features and split the dataset into training and testing sets.
4. Define an MLP classifier using PyTorch and train it with early stopping based on validation AUC.
5. Evaluate the trained model on the test set and save the best model.

Example Parameters are embedded below and can be adjusted as needed.
"""

import os
import sys
import time
import numpy as np
import scanpy as sc
import scvi
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score

# Custom Dataset for Treatment Response Classification
class TreatmentResponseDataset(Dataset):
    def __init__(self, features, labels):
        """
        Args:
            features (np.array): Array of latent features.
            labels (np.array): Array of binary labels (0 or 1).
        """
        self.features = torch.tensor(features, dtype=torch.float32)
        self.labels = torch.tensor(labels, dtype=torch.float32).unsqueeze(1)
    
    def __len__(self):
        return self.features.shape[0]
    
    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]

# Define a Multi-Layer Perceptron (MLP) Classifier
class MLPClassifier(nn.Module):
    def __init__(self, input_dim, hidden_dims=[64, 32], dropout=0.5):
        """
        Args:
            input_dim (int): Dimension of the input latent features.
            hidden_dims (list): List containing sizes of hidden layers.
            dropout (float): Dropout rate to reduce overfitting.
        """
        super(MLPClassifier, self).__init__()
        layers = []
        prev_dim = input_dim
        for hdim in hidden_dims:
            layers.append(nn.Linear(prev_dim, hdim))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
            prev_dim = hdim
        layers.append(nn.Linear(prev_dim, 1))  # Output layer for binary classification
        self.model = nn.Sequential(*layers)
    
    def forward(self, x):
        return self.model(x)

def train_classifier(model, dataloader, criterion, optimizer, device):
    """Train the classifier for one epoch."""
    model.train()
    running_loss = 0.0
    for features, labels in dataloader:
        features, labels = features.to(device), labels.to(device)
        optimizer.zero_grad()
        outputs = model(features)
        loss = criterion(outputs, labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item() * features.size(0)
    epoch_loss = running_loss / len(dataloader.dataset)
    return epoch_loss

def evaluate_classifier(model, dataloader, device):
    """Evaluate the classifier and compute accuracy and ROC AUC."""
    model.eval()
    preds = []
    true_labels = []
    with torch.no_grad():
        for features, labels in dataloader:
            features = features.to(device)
            outputs = model(features)
            preds.append(outputs.cpu().numpy())
            true_labels.append(labels.cpu().numpy())
    preds = np.concatenate(preds, axis=0)
    true_labels = np.concatenate(true_labels, axis=0)
    # Apply sigmoid to convert logits to probabilities
    probs = 1 / (1 + np.exp(-preds))
    binary_preds = (probs >= 0.5).astype(int)
    accuracy = accuracy_score(true_labels, binary_preds)
    try:
        auc = roc_auc_score(true_labels, probs)
    except ValueError:
        auc = 0.0
    return accuracy, auc

def main():
    # Example parameters embedded in the script
    input_h5ad = "output/scanpy/scRNAseq_clustered_treatment.h5ad"  # AnnData file with treatment_response labels
    model_save_path = "output/ml_models/treatment_response_classifier.pt"
    latent_save_path = "output/ml_models/latent_representation.npy"
    label_key = "treatment_response"  # Must be present in adata.obs (values: 0 or 1)
    test_size = 0.2
    random_state = 42
    batch_size = 32
    num_epochs = 50
    learning_rate = 1e-3
    early_stopping_patience = 5

    # Create output directory if it does not exist
    os.makedirs(os.path.dirname(model_save_path), exist_ok=True)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("Using device:", device)

    # Load the AnnData object
    print("Loading AnnData from:", input_h5ad)
    try:
        adata = sc.read_h5ad(input_h5ad)
    except Exception as e:
        print("Error loading AnnData file:", e)
        sys.exit(1)
    
    # Ensure the treatment_response label is present
    if label_key not in adata.obs:
        print(f"Error: '{label_key}' column not found in adata.obs")
        sys.exit(1)
    
    # Setup AnnData for scvi-tools with treatment_response as label
    scvi.data.setup_anndata(adata, labels_key=label_key)
    
    # Train scVI model for latent representation extraction
    print("Training scVI model to obtain latent representation...")
    vae = scvi.model.SCVI(adata, n_latent=10)
    vae.train(max_epochs=100, early_stopping=True)
    
    # Extract latent representation and save for future use
    latent = vae.get_latent_representation()
    print("Latent representation shape:", latent.shape)
    np.save(latent_save_path, latent)
    
    # Extract binary labels from adata.obs
    labels = adata.obs[label_key].astype(int).values
    if latent.shape[0] != len(labels):
        print("Error: The number of latent representations does not match the number of labels.")
        sys.exit(1)
    
    # Split data into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(latent, labels, test_size=test_size, random_state=random_state)
    
    # Create PyTorch datasets and dataloaders
    train_dataset = TreatmentResponseDataset(X_train, y_train)
    test_dataset = TreatmentResponseDataset(X_test, y_test)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    # Define the MLP classifier
    input_dim = latent.shape[1]
    classifier = MLPClassifier(input_dim=input_dim, hidden_dims=[64, 32], dropout=0.5)
    classifier.to(device)
    
    criterion = nn.BCEWithLogitsLoss()
    optimizer = optim.Adam(classifier.parameters(), lr=learning_rate)
    
    best_auc = 0.0
    epochs_no_improve = 0
    best_model_state = None

    print("Starting training of the MLP classifier...")
    for epoch in range(num_epochs):
        train_loss = train_classifier(classifier, train_loader, criterion, optimizer, device)
        train_acc, train_auc = evaluate_classifier(classifier, train_loader, device)
        val_acc, val_auc = evaluate_classifier(classifier, test_loader, device)
        print(f"Epoch {epoch+1}/{num_epochs} - Loss: {train_loss:.4f} - Train Acc: {train_acc:.4f} - Train AUC: {train_auc:.4f} - Val Acc: {val_acc:.4f} - Val AUC: {val_auc:.4f}")
        
        # Early stopping based on validation AUC improvement
        if val_auc > best_auc:
            best_auc = val_auc
            best_model_state = classifier.state_dict()
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= early_stopping_patience:
                print("Early stopping triggered.")
                break
    
    # Load the best model state if available
    if best_model_state is not None:
        classifier.load_state_dict(best_model_state)
    
    # Final evaluation on the test set
    test_acc, test_auc = evaluate_classifier(classifier, test_loader, device)
    print(f"Final Test Accuracy: {test_acc:.4f} - Final Test AUC: {test_auc:.4f}")
    
    # Save the trained classifier model
    torch.save(classifier.state_dict(), model_save_path)
    print("Trained MLP classifier saved to:", model_save_path)

if __name__ == "__main__":
    main()