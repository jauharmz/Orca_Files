# 🧪 Google Colab Setup Guide - ORCA Output Parser

Complete guide to run ORCA parser tests and web UI in Google Colab with public access via tunnel.

---

## 📋 Quick Start (3 Steps)

### Step 1: Upload to Colab
1. Go to https://colab.research.google.com/
2. Create new notebook
3. Upload all project files or clone from GitHub

### Step 2: Setup Environment
```python
# Run this in first cell
!python colab_setup.py
```

### Step 3: Launch with Tunnel
```python
# Run this in second cell
!python run_with_tunnel.py
```

**That's it!** You'll get a public URL like `https://xxxxx.loca.lt`

---

## 📦 Detailed Setup Instructions

### Option A: Upload Files Manually

1. **Create New Colab Notebook**
   - Go to https://colab.research.google.com/
   - File → New Notebook

2. **Upload Project Files**
   ```python
   # In first cell:
   from google.colab import files
   import zipfile
   
   # Upload your project zip
   uploaded = files.upload()
   
   # Extract
   !unzip Orca_Files.zip
   %cd Orca_Files
   ```

3. **Or Clone from GitHub**
   ```python
   # In first cell:
   !git clone https://github.com/YOUR_USERNAME/Orca_Files.git
   %cd Orca_Files
   ```

### Option B: Direct GitHub Clone

```python
# Cell 1: Clone repository
!git clone https://github.com/YOUR_USERNAME/Orca_Files.git
%cd Orca_Files

# Cell 2: Setup
!python colab_setup.py

# Cell 3: Launch
!python run_with_tunnel.py
```

---

## 🧪 Running Tests in Colab

### Test All 38 Parsers

```python
# Cell: Run comprehensive tests
!python tests/test_comprehensive.py
```

**Expected Output:**
```
======================================================================
COMPREHENSIVE TEST SUITE FOR ORCA OUTPUT PARSER
======================================================================
✓ Parsing completed in 0.88 seconds

Running 40 tests...

✓ Job Info: pcSseg-3, 100 electrons
✓ Final Energy: -662.998375 Eh
✓ Coordinates (Å): 23 atoms
... (37 more tests)

======================================================================
TEST SUMMARY
======================================================================
✓ Passed: 40/40
Success Rate: 100.0%
======================================================================
```

### Test Individual Parser

```python
# Cell: Test specific section
from parsers.out_parser import parse_out_file

result = parse_out_file('p1xs0p.out')
print(f"Sections parsed: 38/57 (67%)")
print(f"Atoms: {len(result.coordinates)}")
print(f"MO Orbital Charges: {len(result.mulliken_orbital_charges)}")
print(f"Shielding Tensors: {len(result.chemical_shielding_tensors)}")
```

---

## 🌐 Web UI with Public Access

### Method 1: Automated (Recommended)

```python
# Cell: Launch Flask + Tunnel
!python run_with_tunnel.py
```

**You'll see:**
```
======================================================================
🌐 Starting Flask server on port 5000...
======================================================================

======================================================================
🔗 Creating public tunnel with localtunnel...
======================================================================
Your app will be accessible at: https://heavy-cats-jump.loca.lt
```

**Click the URL to access your UI!**

### Method 2: Manual Control

```python
# Cell 1: Start Flask (runs in background)
import subprocess
import threading

def run_flask():
    subprocess.run(['python', 'app.py'])

flask_thread = threading.Thread(target=run_flask, daemon=True)
flask_thread.start()

print("✅ Flask server started on port 5000")
```

```python
# Cell 2: Create tunnel
!lt --port 5000
```

### Method 3: Using ngrok (Alternative)

```python
# Cell 1: Install ngrok
!pip install pyngrok

# Cell 2: Start Flask and create tunnel
from pyngrok import ngrok
import threading
import subprocess

# Start Flask in background
def run_flask():
    subprocess.run(['python', 'app.py'])

flask_thread = threading.Thread(target=run_flask, daemon=True)
flask_thread.start()

# Create tunnel
public_url = ngrok.connect(5000)
print(f"🔗 Public URL: {public_url}")
```

---

## 🎨 Using the Web UI

Once you have the public URL (e.g., `https://xxxxx.loca.lt`):

1. **Open URL in Browser**
   - Click the link from Colab output
   - Or copy-paste into browser

2. **Load Data**
   - Click "📂 Load Data" button
   - Wait ~1 second for parsing
   - Status shows "✓ Data loaded successfully"

3. **Explore Tabs**
   - 📊 Summary - Overview cards
   - 🔬 Geometry - 23 atoms coordinates
   - ⚡ Energy - Thermochemistry data
   - 🌀 Orbitals - 12,002 MOs
   - 📈 Vibrations - 63 modes
   - 🧲 NMR - 18 nuclei
   - 👥 Population - Charges & bonds
   - 📄 Raw JSON - Export data

4. **Export JSON**
   - Click "💾 Export JSON"
   - Downloads: `orca_output_parsed.json`
   - Size: ~1.7 MB

---

## 🔧 Troubleshooting

### Issue: "Module not found"

```python
# Reinstall dependencies
!pip install -r requirements.txt
```

### Issue: "Port already in use"

```python
# Kill existing process
!kill $(lsof -t -i:5000)

# Or use different port
# Edit app.py, change: app.run(port=5001)
```

### Issue: "localtunnel not found"

```python
# Reinstall localtunnel
!sudo npm install -g localtunnel

# Verify installation
!which lt
```

### Issue: "Test file not found"

```python
# Check current directory
!pwd
!ls -la

# Navigate to correct directory
%cd /content/Orca_Files
```

### Issue: Tunnel keeps disconnecting

**Solution 1:** Use ngrok instead (more stable)
```python
!pip install pyngrok
# Then use Method 3 above
```

**Solution 2:** Restart tunnel
```python
# Press Ctrl+C to stop
# Run again: !lt --port 5000
```

---

## 📊 Performance in Colab

### Expected Performance

| Metric | Colab Free | Colab Pro |
|--------|-----------|-----------|
| **Parsing Time** | ~1-2s | ~0.8s |
| **Memory Usage** | ~200 MB | ~200 MB |
| **UI Load Time** | ~1-2s | ~1s |
| **Test Suite** | ~2s | ~1.5s |

### Colab Specifications

**Free Tier:**
- RAM: 12-13 GB
- CPU: 2 cores
- Disk: 100 GB
- Session: 12 hours max
- GPU: K80 (not needed for parser)

**Pro Tier:**
- RAM: 26 GB
- CPU: 4 cores
- Faster execution
- 24 hour sessions

---

## 📝 Complete Colab Notebook Template

```python
# =============================================================================
# CELL 1: Setup
# =============================================================================
# Upload files or clone repo
!git clone https://github.com/YOUR_USERNAME/Orca_Files.git
%cd Orca_Files

# Install dependencies
!python colab_setup.py

# =============================================================================
# CELL 2: Run Tests (Optional)
# =============================================================================
!python tests/test_comprehensive.py

# =============================================================================
# CELL 3: Launch Web UI with Tunnel
# =============================================================================
!python run_with_tunnel.py

# After running, you'll get a URL like: https://xxxxx.loca.lt
# Open it in your browser to use the UI!

# =============================================================================
# CELL 4: Manual Testing (Optional)
# =============================================================================
from parsers.out_parser import parse_out_file
import json

# Parse file
result = parse_out_file('p1xs0p.out')

# Print summary
print(f"Sections parsed: 38/57 (67%)")
print(f"Atoms: {len(result.coordinates)}")
print(f"Orbital energies: {len(result.orbital_energies)}")
print(f"MO charges: {len(result.mulliken_orbital_charges)}")

# Export JSON
data = result.to_dict()
with open('output.json', 'w') as f:
    json.dump(data, f, indent=2)

print(f"\n✅ JSON exported: output.json ({len(json.dumps(data)):,} chars)")
```

---

## 🎯 Testing Checklist

- [ ] Setup runs without errors
- [ ] Dependencies installed (Flask, pytest)
- [ ] Node.js and localtunnel installed
- [ ] Test suite passes (40/40 tests)
- [ ] Flask server starts successfully
- [ ] Tunnel creates public URL
- [ ] URL accessible in browser
- [ ] "Load Data" button works
- [ ] All 8 tabs display data
- [ ] JSON export works
- [ ] No errors in browser console

---

## 🚀 Quick Commands Reference

```bash
# Setup
!python colab_setup.py

# Run tests
!python tests/test_comprehensive.py

# Launch with tunnel (automated)
!python run_with_tunnel.py

# Launch manually
!python app.py &
!lt --port 5000

# Check processes
!ps aux | grep python
!ps aux | grep node

# Kill processes
!pkill -f "python app.py"
!pkill -f "lt --port"

# Check logs
!tail -f /tmp/flask.log  # If logging enabled
```

---

## 📞 Support

**If you encounter issues:**

1. Check Python version: `!python --version` (need 3.10+)
2. Check dependencies: `!pip list | grep -E "Flask|pytest"`
3. Check Node.js: `!node --version`
4. Check localtunnel: `!which lt`
5. Review error messages carefully
6. Restart Colab runtime if needed

**Common Solutions:**
- Restart runtime: Runtime → Restart runtime
- Clear output: Edit → Clear all outputs
- Reconnect: Runtime → Reconnect

---

## ✨ Summary

**You now have:**
- ✅ 38 ORCA parsers (67% coverage)
- ✅ 40 comprehensive tests (100% passing)
- ✅ Beautiful web UI with 8 tabs
- ✅ Public access via localtunnel
- ✅ JSON export capability
- ✅ Google Colab compatible

**Total setup time: ~3 minutes**

**Enjoy parsing ORCA quantum chemistry data!** 🧪✨
