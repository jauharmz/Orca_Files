"""
Google Colab Setup Script for ORCA Output Parser

This script sets up the ORCA parser and web UI in Google Colab.
It also sets up localtunnel or ngrok for external access.

Usage in Colab:
    !python colab_setup.py
"""

import subprocess
import sys
import os
from pathlib import Path

def run_command(cmd, description=""):
    """Run a command and print output."""
    print(f"\n{'='*70}")
    if description:
        print(f"📦 {description}")
    print(f"{'='*70}")
    print(f"Running: {cmd}\n")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print(result.stderr)
    
    if result.returncode != 0:
        print(f"❌ Command failed with code {result.returncode}")
        return False
    else:
        print(f"✅ {description or 'Command'} completed successfully")
        return True

def main():
    print("\n" + "="*70)
    print("🧪 ORCA OUTPUT PARSER - GOOGLE COLAB SETUP")
    print("="*70)
    
    # Step 1: Install Python dependencies
    run_command(
        "pip install -q Flask==3.0.0 Werkzeug==3.0.1 pytest==7.4.3",
        "Installing Python dependencies"
    )
    
    # Step 2: Install Node.js (for localtunnel)
    print("\n" + "="*70)
    print("📦 Installing Node.js and npm")
    print("="*70)
    
    # Check if node is already installed
    node_check = subprocess.run("which node", shell=True, capture_output=True)
    if node_check.returncode == 0:
        print("✅ Node.js already installed")
    else:
        run_command(
            "curl -fsSL https://deb.nodesource.com/setup_18.x | sudo -E bash - && sudo apt-get install -y nodejs",
            "Installing Node.js 18.x"
        )
    
    # Step 3: Install localtunnel
    run_command(
        "npm install -g localtunnel",
        "Installing localtunnel"
    )
    
    # Step 4: Verify installation
    print("\n" + "="*70)
    print("✅ INSTALLATION VERIFICATION")
    print("="*70)
    
    # Check Python version
    py_version = sys.version.split()[0]
    print(f"Python version: {py_version}")
    
    # Check Flask
    try:
        import flask
        print(f"Flask version: {flask.__version__}")
    except ImportError:
        print("❌ Flask not installed")
    
    # Check Node.js
    node_result = subprocess.run("node --version", shell=True, capture_output=True, text=True)
    if node_result.returncode == 0:
        print(f"Node.js version: {node_result.stdout.strip()}")
    
    # Check npm
    npm_result = subprocess.run("npm --version", shell=True, capture_output=True, text=True)
    if npm_result.returncode == 0:
        print(f"npm version: {npm_result.stdout.strip()}")
    
    # Check lt
    lt_result = subprocess.run("which lt", shell=True, capture_output=True, text=True)
    if lt_result.returncode == 0:
        print(f"localtunnel installed: {lt_result.stdout.strip()}")
    
    print("\n" + "="*70)
    print("🎉 SETUP COMPLETE!")
    print("="*70)
    print("\nNext steps:")
    print("1. Run tests: python tests/test_comprehensive.py")
    print("2. Start Flask app: python app.py")
    print("3. In another terminal: lt --port 5000")
    print("\nOr use the automated script: python run_with_tunnel.py")
    print("="*70 + "\n")

if __name__ == "__main__":
    main()
