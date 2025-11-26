/**
 * Diagnostic Script for ORCA Parser - JavaScript Iteration Error
 *
 * Run this in your browser console to identify which render function is failing.
 *
 * Instructions:
 * 1. Open the ORCA viewer in your browser (http://localhost:5000)
 * 2. Open browser Developer Tools (F12)
 * 3. Go to Console tab
 * 4. Copy and paste this entire script
 * 5. Press Enter
 * 6. The script will identify which function is causing the error
 */

console.log('🔍 ORCA Parser Diagnostic Script');
console.log('================================\n');

// Check if data exists
if (typeof data === 'undefined') {
    console.error('❌ ERROR: `data` variable is not defined. Load data first by clicking "Load Data" button.');
} else {
    console.log('✅ Data variable exists');
    console.log('Data keys:', Object.keys(data).length);

    // Check critical arrays
    const arrayFields = [
        'coordinates',
        'coordinates_au',
        'orbital_energies',
        'frequencies',
        'ir_spectrum',
        'raman_spectrum',
        'scf_energies',
        'mayer_bond_orders',
        'loewdin_bond_orders',
        'scf_iterations'
    ];

    console.log('\n📊 Array Fields Check:');
    arrayFields.forEach(field => {
        const value = data[field];
        const isArray = Array.isArray(value);
        const length = isArray ? value.length : 'N/A';
        const status = isArray ? '✅' : '❌';
        console.log(`${status} ${field}: ${isArray ? 'Array' : typeof value} (length: ${length})`);
    });

    // Check object fields
    const objectFields = [
        'mulliken_charges',
        'loewdin_charges',
        'job_info',
        'dipole_moment',
        'thermochemistry'
    ];

    console.log('\n📦 Object Fields Check:');
    objectFields.forEach(field => {
        const value = data[field];
        const isObject = value && typeof value === 'object' && !Array.isArray(value);
        const keys = isObject ? Object.keys(value).length : 'N/A';
        const status = isObject || value === null ? '✅' : '❌';
        console.log(`${status} ${field}: ${typeof value} (keys: ${keys})`);
    });
}

console.log('\n🧪 Testing Render Functions:');
console.log('================================');

// List of all render functions from renderAllPlots()
const renderFunctions = [
    'renderSCFPlot',
    'renderIRPlot',
    'renderRamanPlot',
    'renderThermoPlot',
    'renderOrbitalPlot',
    'renderNMRPlot',
    'renderCombinedSpectrum',
    'renderUVVisElectric',
    'renderUVVisVelocity',
    'renderTimingBreakdown',
    'renderSCFDetails',
    'renderBasisComposition',
    'renderGridStatistics',
    'renderOptimizationTrajectory',
    'populateVibrationalModes',
    'renderOrbitalChargeHeatmap',
    'renderJCouplingNetwork',
    'renderChemicalShieldingTable',
    'renderHOMOLUMOGapTracker',
    'renderDensityOfStates',
    'renderMullikenOverlapNetwork',
    'renderPolarizabilityVisualization',
    'renderBondOrderComparison',
    'renderChargeComparison',
    'renderEnergyComponents',
    'renderDispersionCorrection',
    'renderOrbitalPopulation',
    'renderInternalCoordinates',
    'renderSCFEfficiency',
    'renderSolvationAnalysis',
    'renderOrbitalDistribution',
    'renderFrequencyAnalysis',
    'renderChargeDistributionPie',
    'renderBondCorrelation',
    'renderEigenvalueSpectrum',
    'renderIRRamanCorrelation',
    'renderMassDistribution',
    'renderDipoleComponents',
    'renderThermoBreakdown',
    'renderAtomStats',
    'renderAngleDistribution',
    'renderElementComposition',
    'renderOptConvergence',
    'renderPopulationCorrelation',
    'renderDistanceAnalysis',
    'renderCoordinationNumbers',
    'renderSCFDashboard',
    'renderEnergyWaterfall',
    'renderPerformanceDashboard',
    'renderMemoryUsage',
    'renderOrbitalGapsLandscape',
    'renderModeClassification',
    'renderNMRDistribution',
    'renderBondStrengthDistribution',
    'renderTimingEfficiency',
    'renderPropertiesSummary',
    'renderMOComposition',
    'renderEnergyDiagram',
    'renderCorrelationMatrix',
    'renderRadialDistribution',
    'renderOptimizationPath3D',
    'renderOverlapMatrix',
    'renderChargeFlow',
    'renderFrequency3D',
    'renderBasisDistribution'
];

let successCount = 0;
let failCount = 0;
let failedFunctions = [];

renderFunctions.forEach(funcName => {
    try {
        if (typeof window[funcName] === 'function') {
            window[funcName]();
            console.log(`✅ ${funcName}`);
            successCount++;
        } else {
            console.warn(`⚠️  ${funcName} - Not found`);
        }
    } catch (error) {
        console.error(`❌ ${funcName} - ERROR: ${error.message}`);
        failCount++;
        failedFunctions.push({
            function: funcName,
            error: error.message,
            stack: error.stack
        });
    }
});

console.log('\n📈 Results:');
console.log(`✅ Successful: ${successCount}`);
console.log(`❌ Failed: ${failCount}`);

if (failedFunctions.length > 0) {
    console.log('\n🐛 Failed Functions Details:');
    failedFunctions.forEach(({function: func, error, stack}) => {
        console.group(`❌ ${func}`);
        console.error('Error:', error);
        console.log('Stack trace:', stack);
        console.groupEnd();
    });

    console.log('\n💡 Solution:');
    console.log('The error is in one of the functions above.');
    console.log('Look for iteration operations (forEach, map, for...of) on non-iterable objects.');
}

console.log('\n================================');
console.log('Diagnostic complete!');
