#!/bin/bash

# Script to install bioinformatics tools on a Google Cloud Vertex AI instance (Debian-based).
# Run this script with sudo privileges, e.g., sudo ./install_bioinfo_tools.sh

echo "Starting bioinformatics tools installation..."
echo "This script will update package lists and install necessary software."
echo "It assumes you are running on a Debian-based Linux system (like Vertex AI default images)."
echo "---------------------------------------------------------------------"

# Exit immediately if a command exits with a non-zero status.
set -e

# 1. Update package lists
echo ""
echo "Updating package lists..."
sudo apt-get update -y
echo "Package lists updated."
echo "---------------------------------------------------------------------"

# 2. Install essential build tools and libraries
# These are often dependencies for compiling or running bioinformatics software.
echo ""
echo "Installing essential build tools and libraries (build-essential, zlib1g-dev, etc.)..."
sudo apt-get install -y \
    build-essential \
    wget \
    curl \
    unzip \
    gzip \
    git \
    zlib1g-dev \
    libncurses5-dev \
    libbz2-dev \
    liblzma-dev \
    autotools-dev \
    autoconf \
    pkg-config
echo "Essential tools and libraries installed."
echo "---------------------------------------------------------------------"

# 3. Install core utilities (most should be present, but good to ensure)
# awk (gawk), grep, cut, sort, uniq are part of coreutils or gawk.
echo ""
echo "Ensuring core utilities are present (gawk, coreutils)..."
sudo apt-get install -y gawk coreutils
echo "Core utilities checked/installed."
echo "---------------------------------------------------------------------"

# 4. Install Samtools
echo ""
echo "Installing Samtools..."
sudo apt-get install -y samtools
echo "Samtools installed."
samtools --version # Verify installation
echo "---------------------------------------------------------------------"

# 5. Install Minimap2
echo ""
echo "Installing Minimap2..."
sudo apt-get install -y minimap2
echo "Minimap2 installed."
minimap2 --version # Verify installation
echo "---------------------------------------------------------------------"

# 6. Install Seqtk
# Seqtk might be in apt repositories, if not, we'll download and compile.
echo ""
echo "Attempting to install Seqtk from apt..."
if sudo apt-get install -y seqtk; then
    echo "Seqtk installed successfully from apt."
else
    echo "Seqtk not found in apt or installation failed. Attempting to install from source..."
    cd /tmp # Go to a temporary directory
    git clone https://github.com/lh3/seqtk.git
    cd seqtk
    make
    sudo mv seqtk /usr/local/bin/
    cd / # Go back to root or home
    rm -rf /tmp/seqtk # Clean up
    echo "Seqtk installed from source to /usr/local/bin/."
fi
seqtk # Verify installation (should show usage)
echo "---------------------------------------------------------------------"

# 7. Install SeqKit
# SeqKit is typically installed by downloading the pre-compiled binary.
echo ""
echo "Installing SeqKit..."
# Get the latest release URL from GitHub (you might want to check for the very latest version)
# For simplicity, using a known recent version. Adjust if needed.
SEQKIT_VERSION="2.10.0" # Example version, check for latest
SEQKIT_ARCH=$(dpkg --print-architecture) # e.g., amd64, arm64

if [ "$SEQKIT_ARCH" == "amd64" ]; then
    SEQKIT_DOWNLOAD_URL="https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_amd64.tar.gz"
elif [ "$SEQKIT_ARCH" == "arm64" ]; then
    SEQKIT_DOWNLOAD_URL="https://github.com/shenwei356/seqkit/releases/download/v${SEQKIT_VERSION}/seqkit_linux_arm64.tar.gz"
else
    echo "Unsupported architecture for SeqKit binary: $SEQKIT_ARCH. Please install manually."
    exit 1
fi

echo "Downloading SeqKit binary from ${SEQKIT_DOWNLOAD_URL}..."
cd /tmp
wget -q "${SEQKIT_DOWNLOAD_URL}" -O seqkit.tar.gz
tar -xzf seqkit.tar.gz
sudo mv seqkit /usr/local/bin/
rm seqkit.tar.gz # Clean up
cd /
echo "SeqKit installed to /usr/local/bin/."
seqkit version # Verify installation
echo "---------------------------------------------------------------------"

echo ""
echo "All specified bioinformatics tools installation process completed."
echo "Please verify each tool by running its --version command or by checking its usage."
echo "Example verifications:"
echo "  samtools --version"
echo "  minimap2 --version"
echo "  seqtk"
echo "  seqkit version"
