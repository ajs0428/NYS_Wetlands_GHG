#!/usr/bin/env python
# coding: utf-8

# In[6]:


import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import json
import time
import sys

# Add script directory to Python path
script_dir = Path("Python_Code_Analysis/DL_Implement/")
sys.path.insert(0, str(script_dir))

# Import modules
from NYS_04_dataset import get_dataloaders, find_patch_files, load_and_merge_metadata
from NYS_05_unet_model import UNet

# === LOAD DATA ===
print("Loading data...")
train_loader, val_loader, metadata = get_dataloaders(
    data_dir, 
    cluster_id=cluster_id, 
    huc_id=huc_id, 
    batch_size=batch_size
)

print(f"\nDataset Summary:")
print(f"  Training batches: {len(train_loader)}")
print(f"  Validation batches: {len(val_loader)}")
print(f"  in_channels: {metadata['in_channels']}")
print(f"  num_classes: {metadata['num_classes']}")
print(f"  band_names: {metadata['band_names']}")
if "hucs_included" in metadata:
    print(f"  HUCs included: {metadata['hucs_included']}")


# In[7]:


def train_one_epoch(model, train_loader, criterion, optimizer, device):
    """Train for one epoch and return average loss."""
    model.train()
    running_loss = 0.0

    for batch_idx, (X, y) in enumerate(train_loader):
        X, y = X.to(device), y.to(device)

        # Forward pass
        optimizer.zero_grad()
        outputs = model(X)
        loss = criterion(outputs, y)

        # Backward pass
        loss.backward()
        optimizer.step()

        running_loss += loss.item()

        # Progress update every 10 batches
        if (batch_idx + 1) % 10 == 0:
            print(f"    Batch {batch_idx + 1}/{len(train_loader)}, Loss: {loss.item():.4f}")

    return running_loss / len(train_loader)


def validate(model, val_loader, criterion, device, metadata):
    """Validate and return loss plus per-class accuracy."""
    model.eval()
    running_loss = 0.0
    num_classes = metadata["num_classes"]
    class_names = metadata["class_names"]

    # Track correct predictions per class
    correct_per_class = torch.zeros(num_classes)
    total_per_class = torch.zeros(num_classes)

    with torch.no_grad():
        for X, y in val_loader:
            X, y = X.to(device), y.to(device)

            outputs = model(X)
            loss = criterion(outputs, y)
            running_loss += loss.item()

            # Get predictions
            preds = torch.argmax(outputs, dim=1)

            # Per-class accuracy
            for c in range(num_classes):
                mask = (y == c)
                total_per_class[c] += mask.sum().item()
                correct_per_class[c] += ((preds == c) & mask).sum().item()

    avg_loss = running_loss / len(val_loader)

    # Calculate per-class accuracy
    class_acc = {}
    for c in range(num_classes):
        if total_per_class[c] > 0:
            class_acc[class_names[c]] = correct_per_class[c] / total_per_class[c]
        else:
            class_acc[class_names[c]] = 0.0

    # Overall accuracy
    overall_acc = correct_per_class.sum() / total_per_class.sum()

    return avg_loss, overall_acc.item(), class_acc


# In[8]:


def compute_class_weights(data_dir, cluster_id, huc_id, class_names):
    """Compute class weights from training data using inverse frequency."""
    files = find_patch_files(data_dir, cluster_id, huc_id)

    print("Loading y_train files for class weight computation...")
    y_train_list = [np.load(f) for f in files['y_train']]
    y_train = np.concatenate(y_train_list, axis=0)
    print(f"  Combined y_train shape: {y_train.shape}")

    # Count pixels per class
    classes, counts = np.unique(y_train, return_counts=True)
    total = counts.sum()

    # Compute inverse frequency weights
    frequencies = counts / total
    weights = 1.0 / frequencies
    weights = weights / weights.min()  # Normalize so smallest weight is 1.0

    print("\nClass distribution and weights:")
    for c, count, w in zip(classes, counts, weights):
        pct = count / total * 100
        print(f"  {class_names[c]}: {count:,} pixels ({pct:.2f}%) -> weight: {w:.2f}")

    return torch.tensor(weights, dtype=torch.float32)


def main():
    # Device selection
    device = torch.device("cuda" if torch.cuda.is_available() else 
                          "mps" if torch.backends.mps.is_available() else "cpu")
    print(f"Using device: {device}")

    # === COMPUTE CLASS WEIGHTS ===
    class_weights = compute_class_weights(
        data_dir, cluster_id, huc_id, metadata["class_names"]
    )

    # === CREATE MODEL ===
    print("\nInitializing model...")
    model = UNet(
        in_channels=metadata["in_channels"],
        num_classes=metadata["num_classes"],
        base_filters=base_filters
    )
    model = model.to(device)

    total_params = sum(p.numel() for p in model.parameters())
    print(f"Total parameters: {total_params:,}")

    # === LOSS AND OPTIMIZER ===
    class_weights = class_weights.to(device)
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    # === TRAINING LOOP ===
    print("\nStarting training...")
    print("=" * 60)

    best_val_loss = float('inf')
    history = {'train_loss': [], 'val_loss': [], 'val_acc': []}

    for epoch in range(num_epochs):
        epoch_start = time.time()
        print(f"\nEpoch {epoch + 1}/{num_epochs}")
        print("-" * 40)

        # Train
        train_loss = train_one_epoch(model, train_loader, criterion, optimizer, device)

        # Validate
        val_loss, val_acc, class_acc = validate(model, val_loader, criterion, device, metadata)

        epoch_time = time.time() - epoch_start

        # Log results
        print(f"\n  Train Loss: {train_loss:.4f}")
        print(f"  Val Loss:   {val_loss:.4f}")
        print(f"  Val Acc:    {val_acc:.4f}")
        print(f"  Time:       {epoch_time:.1f}s")
        print("  Per-class accuracy:")
        for name, acc in class_acc.items():
            print(f"    {name}: {acc:.4f}")

        # Save history
        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['val_acc'].append(val_acc)

        # Save best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'val_loss': val_loss,
                'val_acc': val_acc,
                'metadata': metadata,
                'config': {
                    'cluster_id': cluster_id,
                    'huc_id': huc_id,
                    'base_filters': base_filters,
                    'learning_rate': learning_rate,
                }
            }, output_dir / "best_model.pth")
            print("  [Saved new best model]")

    # === SAVE FINAL MODEL AND HISTORY ===
    torch.save({
        'epoch': num_epochs,
        'model_state_dict': model.state_dict(),
        'optimizer_state_dict': optimizer.state_dict(),
        'val_loss': val_loss,
        'val_acc': val_acc,
        'metadata': metadata,
        'config': {
            'cluster_id': cluster_id,
            'huc_id': huc_id,
            'base_filters': base_filters,
            'learning_rate': learning_rate,
        }
    }, output_dir / "final_model.pth")

    np.save(output_dir / "training_history.npy", history)

    print("\n" + "=" * 60)
    print("Training complete!")
    print(f"Best validation loss: {best_val_loss:.4f}")
    print(f"Models saved to: {output_dir}")

    return history


# Run training (comment out to prevent execution)
history = main()

