"""
Flask Web Application for ORCA Output Viewer

Provides a web interface to visualize parsed ORCA quantum chemistry data.

Usage:
    python app.py
    # Then open http://localhost:5000 in your browser
"""

from flask import Flask, render_template, jsonify, request, send_from_directory
from pathlib import Path
import json
import sys

# Add parsers to path
sys.path.insert(0, str(Path(__file__).parent))
from parsers.out_parser import parse_out_file

app = Flask(__name__)

# Global variable to store parsed data
parsed_data = None
test_file = Path(__file__).parent / "p1xs0p.out"


@app.route('/')
def index():
    """Serve the main UI page."""
    return render_template('index.html')


@app.route('/api/parse', methods=['POST'])
def parse_file():
    """Parse ORCA output file and return JSON data."""
    global parsed_data

    try:
        # For now, parse the test file
        # In future, could upload custom files
        result = parse_out_file(str(test_file))
        parsed_data = result.to_dict()

        return jsonify({
            'success': True,
            'message': f'Parsed {test_file.name} successfully',
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
