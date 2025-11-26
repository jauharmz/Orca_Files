"""
Flask Web Application for ORCA Output Viewer

Provides a web interface to visualize parsed ORCA quantum chemistry data.

Usage:
    python app.py
    # Then open http://localhost:5000 in your browser
"""

from flask import Flask, render_template, jsonify, request, send_from_directory, Response
from pathlib import Path
from werkzeug.utils import secure_filename
import json
import sys
import os
import tempfile
import numpy as np

# Add parsers to path
sys.path.insert(0, str(Path(__file__).parent))
from parsers.out_parser import parse_out_file

# Import Plotly visualization modules
from ir_visualization import create_stacked_ir_plot, create_dual_ir_plot
from raman_visualization import create_stacked_raman_plot, create_single_raman_with_regions
from absorption_visualization import create_absorption_comparison
from emission_visualization import create_emission_spectrum, create_absorption_emission_mirror
from orbital_visualization import create_orbital_comparison

app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100MB max file size
app.config['UPLOAD_FOLDER'] = tempfile.gettempdir()

# Allowed file extensions
ALLOWED_EXTENSIONS = {'out', 'txt'}

# Global variables
parsed_data = None
current_filename = None
test_file = Path(__file__).parent / "p1xs0p.out"

def allowed_file(filename):
    """Check if file extension is allowed."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/')
def index():
    """Serve the main UI page."""
    return render_template('index.html')


@app.route('/api/upload', methods=['POST'])
def upload_file():
    """Upload and parse ORCA output file."""
    global parsed_data, current_filename

    # Check if file was uploaded
    if 'file' not in request.files:
        return jsonify({
            'success': False,
            'message': 'No file uploaded'
        }), 400

    file = request.files['file']

    # Check if filename is empty
    if file.filename == '':
        return jsonify({
            'success': False,
            'message': 'No file selected'
        }), 400

    # Check if file type is allowed
    if not allowed_file(file.filename):
        return jsonify({
            'success': False,
            'message': f'File type not allowed. Allowed types: {", ".join(ALLOWED_EXTENSIONS)}'
        }), 400

    try:
        # Save file temporarily
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)

        # Parse the file
        print(f"\n📄 Parsing uploaded file: {filename}")
        result = parse_out_file(filepath)
        parsed_data = result.to_dict()
        current_filename = filename

        # Clean up temp file
        os.remove(filepath)

        return jsonify({
            'success': True,
            'message': f'Successfully parsed {filename}',
            'filename': filename,
            'data': parsed_data
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'message': f'Error parsing file: {str(e)}'
        }), 500


@app.route('/api/parse', methods=['POST'])
def parse_default_file():
    """Parse the default test file."""
    global parsed_data, current_filename

    try:
        result = parse_out_file(str(test_file))
        parsed_data = result.to_dict()
        current_filename = test_file.name

        return jsonify({
            'success': True,
            'message': f'Parsed {test_file.name} successfully',
            'filename': test_file.name,
            'data': parsed_data
        })

    except Exception as e:
        return jsonify({
            'success': False,
            'message': f'Error parsing file: {str(e)}'
        }), 500


@app.route('/api/data')
def get_data():
    """Return cached parsed data."""
    if parsed_data is None:
        return jsonify({
            'success': False,
            'message': 'No data loaded. Click "Load Data" first.'
        }), 404

    return jsonify({
        'success': True,
        'data': parsed_data
    })


@app.route('/api/summary')
def get_summary():
    """Return quick summary statistics."""
    if parsed_data is None:
        return jsonify({'success': False, 'message': 'No data loaded'}), 404

    summary = {
        'file_info': {
            'basis_set': parsed_data['job_info']['basis_set'],
            'charge': parsed_data['job_info']['charge'],
            'multiplicity': parsed_data['job_info']['multiplicity'],
            'num_electrons': parsed_data['job_info']['num_electrons'],
        },
        'energy': {
            'final_energy': parsed_data.get('final_energy'),
            'gibbs_free_energy': parsed_data.get('thermochemistry', {}).get('gibbs_free_energy'),
        },
        'geometry': {
            'num_atoms': len(parsed_data.get('coordinates', [])),
            'dipole_moment': parsed_data.get('dipole_moment', {}).get('magnitude_debye'),
        },
        'spectroscopy': {
            'num_frequencies': len(parsed_data.get('frequencies', [])),
            'num_raman_modes': len(parsed_data.get('raman_spectrum', [])),
            'num_nmr_shifts': len(parsed_data.get('nmr_data', {}).get('chemical_shifts', [])),
        },
        'population': {
            'num_mulliken_charges': len(parsed_data.get('mulliken_charges', {})),
            'num_bond_orders': len(parsed_data.get('mayer_bond_orders', [])),
        },
        'sections_parsed': 34,
        'total_sections': 57,
        'coverage_percent': 60
    }

    return jsonify({
        'success': True,
        'summary': summary
    })


@app.route('/api/export/json')
def export_json():
    """Export parsed data as JSON file."""
    if parsed_data is None:
        return jsonify({'success': False, 'message': 'No data loaded'}), 404

    filename = current_filename.replace('.out', '_parsed.json') if current_filename else 'orca_data.json'

    return Response(
        json.dumps(parsed_data, indent=2),
        mimetype='application/json',
        headers={'Content-Disposition': f'attachment; filename={filename}'}
    )


@app.route('/api/visualize/ir')
def visualize_ir():
    """Generate interactive IR spectrum visualization."""
    if parsed_data is None:
        return jsonify({'success': False, 'message': 'No data loaded'}), 404

    try:
        frequencies = parsed_data.get('frequencies', [])
        ir_intensities = [mode.get('ir_intensity', 0) for mode in parsed_data.get('vibrational_modes', [])]

        if not frequencies or not ir_intensities:
            return jsonify({'success': False, 'message': 'No IR data available'}), 404

        # Create IR visualization
        fig = create_dual_ir_plot(
            np.array(frequencies),
            np.array(ir_intensities),
            title=f"IR Spectrum - {current_filename}"
        )

        # Return as HTML
        html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        return Response(html, mimetype='text/html')

    except Exception as e:
        return jsonify({'success': False, 'message': f'Error creating visualization: {str(e)}'}), 500


@app.route('/api/visualize/raman')
def visualize_raman():
    """Generate interactive Raman spectrum visualization."""
    if parsed_data is None:
        return jsonify({'success': False, 'message': 'No data loaded'}), 404

    try:
        raman_data = parsed_data.get('raman_spectrum', [])

        if not raman_data:
            return jsonify({'success': False, 'message': 'No Raman data available'}), 404

        frequencies = np.array([mode['frequency'] for mode in raman_data])
        activities = np.array([mode['raman_activity'] for mode in raman_data])

        # Create Raman visualization
        fig = create_single_raman_with_regions(
            frequencies,
            activities,
            title=f"Raman Spectrum - {current_filename}"
        )

        # Return as HTML
        html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        return Response(html, mimetype='text/html')

    except Exception as e:
        return jsonify({'success': False, 'message': f'Error creating visualization: {str(e)}'}), 500


@app.route('/api/visualize/orbitals')
def visualize_orbitals():
    """Generate interactive orbital energy diagram."""
    if parsed_data is None:
        return jsonify({'success': False, 'message': 'No data loaded'}), 404

    try:
        orbitals = parsed_data.get('orbital_energies', [])

        if not orbitals:
            return jsonify({'success': False, 'message': 'No orbital data available'}), 404

        # Create orbital visualization
        dataset = {'orbital_energies': orbitals}
        fig = create_orbital_comparison(
            [dataset],
            [current_filename],
            title=f"Orbital Energies - {current_filename}"
        )

        # Return as HTML
        html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        return Response(html, mimetype='text/html')

    except Exception as e:
        return jsonify({'success': False, 'message': f'Error creating visualization: {str(e)}'}), 500


@app.route('/api/visualize/absorption')
def visualize_absorption():
    """Generate interactive absorption spectrum visualization."""
    if parsed_data is None:
        return jsonify({'success': False, 'message': 'No data loaded'}), 404

    try:
        transitions = parsed_data.get('tddft_transitions', [])

        if not transitions:
            return jsonify({'success': False, 'message': 'No absorption data available'}), 404

        # Create absorption visualization
        dataset = {'transitions': transitions}
        fig = create_absorption_comparison(
            [dataset],
            [current_filename],
            title=f"UV-Vis Absorption Spectrum - {current_filename}"
        )

        # Return as HTML
        html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        return Response(html, mimetype='text/html')

    except Exception as e:
        return jsonify({'success': False, 'message': f'Error creating visualization: {str(e)}'}), 500


if __name__ == '__main__':
    print("=" * 70)
    print("ORCA Output Viewer - Web Interface")
    print("=" * 70)
    print(f"\nTest file: {test_file}")
    print(f"File exists: {test_file.exists()}")
    print("\nStarting Flask server...")
    print("Open your browser to: http://localhost:5000")
    print("=" * 70)

    app.run(debug=True, host='0.0.0.0', port=5000)
