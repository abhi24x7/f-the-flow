/**
 * Friction Factor Correlations for Pipe Flow
 * All correlations return the Darcy friction factor (f)
 * 
 * Parameters:
 * - eOverD: Relative roughness (e/D)
 * - re: Reynolds number
 * 
 * Returns: Friction factor (f) or null if invalid
 */

const frictionFactorCorrelations = {
    
    /**
     * Swamee-Jain Correlation (1976)
     * Valid for: 5000 < Re < 10^8, 10^-6 < e/D < 10^-2
     * Accuracy: ±1%
     */
    swameeJain: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const term = Math.log10(eOverD / 3.7 + 5.74 / Math.pow(re, 0.9));
        return 0.25 / Math.pow(term, 2);
    },
    
    /**
     * Colebrook Equation (1939)
     * The implicit equation solved iteratively
     * Most accurate but computationally intensive
     */
    colebrook: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        let f = 0.02; // Initial guess
        let iterations = 0;
        const maxIterations = 100;
        const tolerance = 1e-10;
        
        while (iterations < maxIterations) {
            const sqrtF = Math.sqrt(f);
            const term = eOverD / 3.7 + 2.51 / (re * sqrtF);
            const fNew = 0.25 / Math.pow(Math.log10(term), 2);
            
            if (Math.abs(fNew - f) < tolerance) {
                return fNew;
            }
            
            f = fNew;
            iterations++;
        }
        
        return f; // Return last iteration if not converged
    },
    
    /**
     * Haaland Correlation (1983)
     * Explicit approximation to Colebrook
     * Valid for: 4000 < Re < 10^8, 10^-6 < e/D < 5×10^-2
     */
    haaland: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const term1 = Math.pow(eOverD / 3.7, 1.11);
        const term2 = 6.9 / re;
        const denominator = Math.log10(term1 + term2);
        return Math.pow(-1.8 * denominator, -2);
    },
    
    /**
     * Chen Correlation (1979)
     * Explicit approximation with high accuracy
     * Valid for: 4000 < Re < 4×10^8, 10^-7 < e/D < 5×10^-2
     */
    chen: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const A = Math.log10(
            eOverD / 3.7065 - 
            5.0452 / re * Math.log10(
                Math.pow(eOverD / 2.8257, 1.1098) + 
                Math.pow(5.8506 / re, 0.8981)
            )
        );
        return Math.pow(-2 * A, -2);
    },
    
    /**
     * Zigrang-Sylvester Correlation (1982)
     * Two-step explicit approximation
     * Valid for: 4000 < Re < 10^8, 4×10^-5 < e/D < 5×10^-2
     */
    zigrangSylvester: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const A = -0.8686 * Math.log10(eOverD / 3.7 + 5.74 / Math.pow(re, 0.9));
        const B = A - 0.8686 * Math.log10(eOverD / 3.7 + 2.51 * A / re);
        return Math.pow(A - (Math.pow(A, 2) - B) / (2 * A - B), -2);
    },
    
    /**
     * Goudar-Sonnad Correlation (2008)
     * Most accurate explicit approximation
     * Valid for: 4000 < Re < 10^8, 10^-6 < e/D < 5×10^-2
     */
    goudarSonnad: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const a = 2 / Math.log(10);
        const b = eOverD / 3.7;
        const d = Math.log(10) / 5.02 * re;
        const s = b * d + Math.log(d);
        const q = Math.pow(s, Math.pow(s, 1) / (Math.pow(s, 1) + 1));
        const g = b * d + Math.log(d / q);
        const z = Math.log(q / g);
        const dLA = g / (g + 1) * z;
        const dCF = dLA * (1 + z / 2 / Math.pow(g + 1, 2));
        return a * Math.pow(Math.log(d / q) + dCF, -2);
    },
    
    /**
     * Serghides Correlation (1984)
     * Three-step explicit method
     * High accuracy with reasonable computational cost
     */
    serghides: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const A = -2 * Math.log10(eOverD / 3.7 + 12 / re);
        const B = -2 * Math.log10(eOverD / 3.7 + 2.51 * A / re);
        const C = -2 * Math.log10(eOverD / 3.7 + 2.51 * B / re);
        
        return Math.pow(A - Math.pow(B - A, 2) / (C - 2 * B + A), -2);
    },
    
    /**
     * Wood Correlation (1966)
     * Simple explicit approximation
     * Less accurate but very fast
     */
    wood: function(eOverD, re) {
        if (re < 4000) return 64 / re; // Laminar flow
        
        const a = 0.53 * eOverD + 0.094 * Math.pow(eOverD, 0.225) + 88 * Math.pow(eOverD, 0.44);
        const b = Math.pow(re, 1.62 * Math.pow(eOverD, 0.134));
        return 0.094 * Math.pow(eOverD, 0.225) + 0.53 * eOverD + 88 * Math.pow(eOverD, 0.44) * Math.pow(re, -1.62 * Math.pow(eOverD, 0.134));
    }
};

/**
 * Utility function to get all available correlation names
 * @returns {Array} Array of correlation names
 */
function getAvailableCorrelations() {
    return Object.keys(frictionFactorCorrelations);
}

/**
 * Utility function to calculate friction factor with validation
 * @param {string} correlationName - Name of the correlation to use
 * @param {number} eOverD - Relative roughness
 * @param {number} re - Reynolds number
 * @returns {Object} Result object with value, error, and info
 */
function calculateFrictionFactor(correlationName, eOverD, re) {
    // Input validation
    if (!frictionFactorCorrelations[correlationName]) {
        return {
            value: null,
            error: `Correlation '${correlationName}' not found`,
            info: null
        };
    }
    
    if (isNaN(eOverD) || isNaN(re)) {
        return {
            value: null,
            error: 'Invalid input: e/D and Re must be numbers',
            info: null
        };
    }
    
    if (eOverD < 0) {
        return {
            value: null,
            error: 'e/D must be non-negative',
            info: null
        };
    }
    
    if (re <= 0) {
        return {
            value: null,
            error: 'Reynolds number must be positive',
            info: null
        };
    }
    
    try {
        const f = frictionFactorCorrelations[correlationName](eOverD, re);
        
        if (f === null || isNaN(f) || f <= 0) {
            return {
                value: null,
                error: 'Calculation failed - check input ranges',
                info: null
            };
        }
        
        // Add flow regime info
        let flowInfo = '';
        if (re < 2300) {
            flowInfo = 'Laminar flow (Re < 2300)';
        } else if (re < 4000) {
            flowInfo = 'Transitional flow (2300 < Re < 4000)';
        } else {
            flowInfo = 'Turbulent flow (Re > 4000)';
        }
        
        return {
            value: f,
            error: null,
            info: flowInfo
        };
        
    } catch (error) {
        return {
            value: null,
            error: 'Calculation error: ' + error.message,
            info: null
        };
    }
}

// Export for use in other files (if using modules)
// export { frictionFactorCorrelations, calculateFrictionFactor, getAvailableCorrelations };