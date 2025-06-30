#!/bin/bash

echo "🔧 Installing requirements..."
pip install biopython tqdm matplotlib

echo "📁 Creating genome_data folder..."
mkdir -p genome_data

echo "✅ All set! Now run: python main.py"
