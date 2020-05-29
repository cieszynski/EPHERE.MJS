/* EPHERE.MJS */

const log = console.log;

const L_TERMS = [
    [
        [175347046.0, 0, 0],
        [3341656.0, 4.6692568, 6283.07585],
        [34894.0, 4.6261, 12566.1517],
        [3497.0, 2.7441, 5753.3849],
        [3418.0, 2.8289, 3.5231],
        [3136.0, 3.6277, 77713.7715],
        [2676.0, 4.4181, 7860.4194],
        [2343.0, 6.1352, 3930.2097],
        [1324.0, 0.7425, 11506.7698],
        [1273.0, 2.0371, 529.691],
        [1199.0, 1.1096, 1577.3435],
        [990, 5.233, 5884.927],
        [902, 2.045, 26.298],
        [857, 3.508, 398.149],
        [780, 1.179, 5223.694],
        [753, 2.533, 5507.553],
        [505, 4.583, 18849.228],
        [492, 4.205, 775.523],
        [357, 2.92, 0.067],
        [317, 5.849, 11790.629],
        [284, 1.899, 796.298],
        [271, 0.315, 10977.079],
        [243, 0.345, 5486.778],
        [206, 4.806, 2544.314],
        [205, 1.869, 5573.143],
        [202, 2.458, 6069.777],
        [156, 0.833, 213.299],
        [132, 3.411, 2942.463],
        [126, 1.083, 20.775],
        [115, 0.645, 0.98],
        [103, 0.636, 4694.003],
        [102, 0.976, 15720.839],
        [102, 4.267, 7.114],
        [99, 6.21, 2146.17],
        [98, 0.68, 155.42],
        [86, 5.98, 161000.69],
        [85, 1.3, 6275.96],
        [85, 3.67, 71430.7],
        [80, 1.81, 17260.15],
        [79, 3.04, 12036.46],
        [75, 1.76, 5088.63],
        [74, 3.5, 3154.69],
        [74, 4.68, 801.82],
        [70, 0.83, 9437.76],
        [62, 3.98, 8827.39],
        [61, 1.82, 7084.9],
        [57, 2.78, 6286.6],
        [56, 4.39, 14143.5],
        [56, 3.47, 6279.55],
        [52, 0.19, 12139.55],
        [52, 1.33, 1748.02],
        [51, 0.28, 5856.48],
        [49, 0.49, 1194.45],
        [41, 5.37, 8429.24],
        [41, 2.4, 19651.05],
        [39, 6.17, 10447.39],
        [37, 6.04, 10213.29],
        [37, 2.57, 1059.38],
        [36, 1.71, 2352.87],
        [36, 1.78, 6812.77],
        [33, 0.59, 17789.85],
        [30, 0.44, 83996.85],
        [30, 2.74, 1349.87],
        [25, 3.16, 4690.48]
    ],
    [
        [628331966747.0, 0, 0],
        [206059.0, 2.678235, 6283.07585],
        [4303.0, 2.6351, 12566.1517],
        [425.0, 1.59, 3.523],
        [119.0, 5.796, 26.298],
        [109.0, 2.966, 1577.344],
        [93, 2.59, 18849.23],
        [72, 1.14, 529.69],
        [68, 1.87, 398.15],
        [67, 4.41, 5507.55],
        [59, 2.89, 5223.69],
        [56, 2.17, 155.42],
        [45, 0.4, 796.3],
        [36, 0.47, 775.52],
        [29, 2.65, 7.11],
        [21, 5.34, 0.98],
        [19, 1.85, 5486.78],
        [19, 4.97, 213.3],
        [17, 2.99, 6275.96],
        [16, 0.03, 2544.31],
        [16, 1.43, 2146.17],
        [15, 1.21, 10977.08],
        [12, 2.83, 1748.02],
        [12, 3.26, 5088.63],
        [12, 5.27, 1194.45],
        [12, 2.08, 4694],
        [11, 0.77, 553.57],
        [10, 1.3, 6286.6],
        [10, 4.24, 1349.87],
        [9, 2.7, 242.73],
        [9, 5.64, 951.72],
        [8, 5.3, 2352.87],
        [6, 2.65, 9437.76],
        [6, 4.67, 4690.48]
    ],
    [
        [52919.0, 0, 0],
        [8720.0, 1.0721, 6283.0758],
        [309.0, 0.867, 12566.152],
        [27, 0.05, 3.52],
        [16, 5.19, 26.3],
        [16, 3.68, 155.42],
        [10, 0.76, 18849.23],
        [9, 2.06, 77713.77],
        [7, 0.83, 775.52],
        [5, 4.66, 1577.34],
        [4, 1.03, 7.11],
        [4, 3.44, 5573.14],
        [3, 5.14, 796.3],
        [3, 6.05, 5507.55],
        [3, 1.19, 242.73],
        [3, 6.12, 529.69],
        [3, 0.31, 398.15],
        [3, 2.28, 553.57],
        [2, 4.38, 5223.69],
        [2, 3.75, 0.98]
    ],
    [
        [289.0, 5.844, 6283.076],
        [35, 0, 0],
        [17, 5.49, 12566.15],
        [3, 5.2, 155.42],
        [1, 4.72, 3.52],
        [1, 5.3, 18849.23],
        [1, 5.97, 242.73]
    ],
    [
        [114.0, 3.142, 0],
        [8, 4.13, 6283.08],
        [1, 3.84, 12566.15]
    ],
    [
        [1, 3.14, 0]
    ]
];
const B_TERMS = [
    [
        [280.0, 3.199, 84334.662],
        [102.0, 5.422, 5507.553],
        [80, 3.88, 5223.69],
        [44, 3.7, 2352.87],
        [32, 4, 1577.34]
    ],
    [
        [9, 3.9, 5507.55],
        [6, 1.73, 5223.69]
    ]
];
const R_TERMS = [
    [
        [100013989.0, 0, 0],
        [1670700.0, 3.0984635, 6283.07585],
        [13956.0, 3.05525, 12566.1517],
        [3084.0, 5.1985, 77713.7715],
        [1628.0, 1.1739, 5753.3849],
        [1576.0, 2.8469, 7860.4194],
        [925.0, 5.453, 11506.77],
        [542.0, 4.564, 3930.21],
        [472.0, 3.661, 5884.927],
        [346.0, 0.964, 5507.553],
        [329.0, 5.9, 5223.694],
        [307.0, 0.299, 5573.143],
        [243.0, 4.273, 11790.629],
        [212.0, 5.847, 1577.344],
        [186.0, 5.022, 10977.079],
        [175.0, 3.012, 18849.228],
        [110.0, 5.055, 5486.778],
        [98, 0.89, 6069.78],
        [86, 5.69, 15720.84],
        [86, 1.27, 161000.69],
        [65, 0.27, 17260.15],
        [63, 0.92, 529.69],
        [57, 2.01, 83996.85],
        [56, 5.24, 71430.7],
        [49, 3.25, 2544.31],
        [47, 2.58, 775.52],
        [45, 5.54, 9437.76],
        [43, 6.01, 6275.96],
        [39, 5.36, 4694],
        [38, 2.39, 8827.39],
        [37, 0.83, 19651.05],
        [37, 4.9, 12139.55],
        [36, 1.67, 12036.46],
        [35, 1.84, 2942.46],
        [33, 0.24, 7084.9],
        [32, 0.18, 5088.63],
        [32, 1.78, 398.15],
        [28, 1.21, 6286.6],
        [28, 1.9, 6279.55],
        [26, 4.59, 10447.39]
    ],
    [
        [103019.0, 1.10749, 6283.07585],
        [1721.0, 1.0644, 12566.1517],
        [702.0, 3.142, 0],
        [32, 1.02, 18849.23],
        [31, 2.84, 5507.55],
        [25, 1.32, 5223.69],
        [18, 1.42, 1577.34],
        [10, 5.91, 10977.08],
        [9, 1.42, 6275.96],
        [9, 0.27, 5486.78]
    ],
    [
        [4359.0, 5.7846, 6283.0758],
        [124.0, 5.579, 12566.152],
        [12, 3.14, 0],
        [9, 3.63, 77713.77],
        [6, 1.87, 5573.14],
        [3, 5.47, 18849.23]
    ],
    [
        [145.0, 4.273, 6283.076],
        [7, 3.92, 12566.15]
    ],
    [
        [4, 2.56, 6283.08]
    ]
];
const PE_TERMS = [
    [-171996, -174.2, 92025, 8.9],
    [-13187, -1.6, 5736, -3.1],
    [-2274, -0.2, 977, -0.5],
    [2062, 0.2, -895, 0.5],
    [1426, -3.4, 54, -0.1],
    [712, 0.1, -7, 0],
    [-517, 1.2, 224, -0.6],
    [-386, -0.4, 200, 0],
    [-301, 0, 129, -0.1],
    [217, -0.5, -95, 0.3],
    [-158, 0, 0, 0],
    [129, 0.1, -70, 0],
    [123, 0, -53, 0],
    [63, 0, 0, 0],
    [63, 0.1, -33, 0],
    [-59, 0, 26, 0],
    [-58, -0.1, 32, 0],
    [-51, 0, 27, 0],
    [48, 0, 0, 0],
    [46, 0, -24, 0],
    [-38, 0, 16, 0],
    [-31, 0, 13, 0],
    [29, 0, 0, 0],
    [29, 0, -12, 0],
    [26, 0, 0, 0],
    [-22, 0, 0, 0],
    [21, 0, -10, 0],
    [17, -0.1, 0, 0],
    [16, 0, -8, 0],
    [-16, 0.1, 7, 0],
    [-15, 0, 9, 0],
    [-13, 0, 7, 0],
    [-12, 0, 6, 0],
    [11, 0, 0, 0],
    [-10, 0, 5, 0],
    [-8, 0, 3, 0],
    [7, 0, -3, 0],
    [-7, 0, 0, 0],
    [-7, 0, 3, 0],
    [-7, 0, 3, 0],
    [6, 0, 0, 0],
    [6, 0, -3, 0],
    [6, 0, -3, 0],
    [-6, 0, 3, 0],
    [-6, 0, 3, 0],
    [5, 0, 0, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [-5, 0, 3, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [-4, 0, 0, 0],
    [3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
    [-3, 0, 0, 0],
];
const Y_TERMS = [
    [0, 0, 0, 0, 1],
    [-2, 0, 0, 2, 2],
    [0, 0, 0, 2, 2],
    [0, 0, 0, 0, 2],
    [0, 1, 0, 0, 0],
    [0, 0, 1, 0, 0],
    [-2, 1, 0, 2, 2],
    [0, 0, 0, 2, 1],
    [0, 0, 1, 2, 2],
    [-2, -1, 0, 2, 2],
    [-2, 0, 1, 0, 0],
    [-2, 0, 0, 2, 1],
    [0, 0, -1, 2, 2],
    [2, 0, 0, 0, 0],
    [0, 0, 1, 0, 1],
    [2, 0, -1, 2, 2],
    [0, 0, -1, 0, 1],
    [0, 0, 1, 2, 1],
    [-2, 0, 2, 0, 0],
    [0, 0, -2, 2, 1],
    [2, 0, 0, 2, 2],
    [0, 0, 2, 2, 2],
    [0, 0, 2, 0, 0],
    [-2, 0, 1, 2, 2],
    [0, 0, 0, 2, 0],
    [-2, 0, 0, 2, 0],
    [0, 0, -1, 2, 1],
    [0, 2, 0, 0, 0],
    [2, 0, -1, 0, 1],
    [-2, 2, 0, 2, 2],
    [0, 1, 0, 0, 1],
    [-2, 0, 1, 0, 1],
    [0, -1, 0, 0, 1],
    [0, 0, 2, -2, 0],
    [2, 0, -1, 2, 1],
    [2, 0, 1, 2, 2],
    [0, 1, 0, 2, 2],
    [-2, 1, 1, 0, 0],
    [0, -1, 0, 2, 2],
    [2, 0, 0, 2, 1],
    [2, 0, 1, 0, 0],
    [-2, 0, 2, 2, 2],
    [-2, 0, 1, 2, 1],
    [2, 0, -2, 0, 1],
    [2, 0, 0, 0, 1],
    [0, -1, 1, 0, 0],
    [-2, -1, 0, 2, 1],
    [-2, 0, 0, 0, 1],
    [0, 0, 2, 2, 1],
    [-2, 0, 2, 0, 1],
    [-2, 1, 0, 2, 1],
    [0, 0, 1, -2, 0],
    [-1, 0, 1, 0, 0],
    [-2, 1, 0, 0, 0],
    [1, 0, 0, 0, 0],
    [0, 0, 1, 2, 0],
    [0, 0, -2, 2, 2],
    [-1, -1, 1, 0, 0],
    [0, 1, 1, 0, 0],
    [0, -1, 1, 2, 2],
    [2, -1, -1, 2, 2],
    [0, 0, 3, 2, 2],
    [2, -1, 0, 2, 2],
];


const PI = Math.PI,
    J2000 = 2451545.0,
    SUN_RADIUS = 0.26667,
    floor = Math.floor,
    round = Math.round,
    pow = Math.pow,
    trunc = Math.trunc,
    sin = Math.sin,
    asin = Math.asin,
    cos = Math.cos,
    acos = Math.acos,
    tan = Math.tan,
    atan = Math.atan,
    atan2 = Math.atan2,
    abs = Math.abs,
    deg = 180 / PI,
    rad = PI / 180;

const frac = (x) => { return x % 1; }

const range = (x, max = 360, min = 0) => {
    return ((max - min) * frac(x / (max - min)))
        + ((max - min) * frac(x / (max - min)) < min) * (max - min)
        - ((max - min) * frac(x / (max - min)) > max) * (max - min);
}

const hhffff = (hh, mm, ss) => { return hh + mm / 60 + ss / 3600; }
const hhmmss = (hhdec) => { return [trunc(hhdec), floor(hhdec % 1 / (1 / 60)), (hhdec % 1 / (1 / 60)) % 1 / (1 / 60)]; }

// https://de.wikipedia.org/wiki/Julianisches_Datum#Zeitma%C3%9Fe_in_Software
export const JD = (y, m, d, hh = 0, mm = 0, ss = 0) => { return 2440587.5 + Date.UTC(y, m - 1, d, hh, mm, ss) / 86400000; }

// The Modified Julian Date (MJD)
export const MJD = (y, m, d, hh = 0, mm = 0, ss = 0) => { return JD(y, m, d, hh, mm, ss) - 2400000.5; }

// Julian century (T or JC)
const T = (y, m, d, hh = 0, mm = 0, ss = 0) => { return (JD(y, m, d, hh, mm, ss) - J2000) / 36525; }, JC = T;
const JDE = (y, m, d, hh, mm, ss) => { return JD(y, m, d, hh, mm, ss) + (/* deltaT(y)  */64.8 / 86400); } // TODO
const JCE = (y, m, d, hh, mm, ss) => { return (JDE(y, m, d, hh, mm, ss) - J2000) / 36525; }

// Julian Ephemeris Millennium (JME) for the 2000 standard epoch
const JME = (y, m, d, hh, mm, ss) => { return JCE(y, m, d, hh, mm, ss) / 10; }
// https://books.google.de/books?id=TppADwAAQBAJ&pg=PA361&lpg=PA361&dq=%22Julian+Ephemeris+Day%22&source=bl&ots=B02n6VIiMq&sig=ACfU3U3dCUFWZO3iTAflooIUS0wdbgbJRQ&hl=de&sa=X&ved=2ahUKEwjRg5v8zdLpAhXQzqQKHfX9B9wQ6AEwBnoECAgQAQ#v=onepage&q=%22Julian%20Ephemeris%20Day%22&f=false
const deltaT = (y) => {
    const tau = (y - 2000);
    if (y >= 1986 && y <= 2005) {
        return 63.86
            + 0.3345 * tau
            - 0.060374 * pow(tau, 2)
            + 0.0017275 * pow(tau, 3)
            + 0.000651814 * pow(tau, 4)
            + 0.00002373599 * pow(tau, 5);
    }

    if (y > 2005 && y <= 2050) {
        return 62.92
            + 0.32217 * tau
            + 0.005589 * pow(tau, 2);
    }

    if (y > 2050 && y < 2150) {
        return - 20 + 32 * pow(round((y - 1820) / 100), 2) - 0.5629 * (2150 - y);
    }

    throw Error(`deltaT(y: ${y}) must between 1986 & 2150!`)
}

// Earth heliocentric longitude [degrees]
const earth_heliocentric_longitude = (jme) => {
    let result = 0;
    L_TERMS.forEach((value, index) => {
        result += (value.map(([A, B, C]) => A * Math.cos(B + C * jme)).reduce((x, y) => x + y)) * pow(jme, index);
    })
    return range(result / 1.0e8 * deg, 360);
}

// Earth heliocentric latitude [degrees]
const earth_heliocentric_latitude = (jme) => {
    let result = 0;
    B_TERMS.forEach((value, index) => {
        result += (value.map(([A, B, C]) => A * Math.cos(B + C * jme)).reduce((x, y) => x + y)) * pow(jme, index);
    })
    return result / 1.0e8 * deg;
}

// Earth radius vector [Astronomical Units, AU
const earth_radius_vector = (jme) => {
    let result = 0;
    R_TERMS.forEach((value, index) => {
        result += (value.map(([A, B, C]) => A * Math.cos(B + C * jme)).reduce((x, y) => x + y)) * pow(jme, index);
    })
    return result / 1.0e8;
}

// Geocentric longitude [degrees]
const geocentric_longitude = (L) => {
    let Theta = L + 180.0;
    if (Theta >= 360.0) Theta -= 360.0;
    return Theta;
}
// Geocentric latitude [degrees]
const geocentric_latitude = (B) => { return -B; }

const third_order_polynomial = (a, b, c, d, x) => { return ((a * x + b) * x + c) * x + d; }
const mean_elongation_moon_sun = (jce) => third_order_polynomial(1.0 / 189474.0, -0.0019142, 445267.11148, 297.85036, jce);
const mean_anomaly_sun = (jce) => third_order_polynomial(-1.0 / 300000.0, -0.0001603, 35999.05034, 357.52772, jce);
const mean_anomaly_moon = (jce) => third_order_polynomial(1.0 / 56250.0, 0.0086972, 477198.867398, 134.96298, jce);
const argument_latitude_moon = (jce) => third_order_polynomial(1.0 / 327270.0, -0.0036825, 483202.017538, 93.27191, jce);
const ascending_longitude_moon = (jce) => third_order_polynomial(1.0 / 450000.0, 0.0020708, -1934.136261, 125.04452, jce);

const nutation_longitude_and_obliquity = (jce) => {
    let data = [
        mean_elongation_moon_sun(jce),
        mean_anomaly_sun(jce),
        mean_anomaly_moon(jce),
        argument_latitude_moon(jce),
        ascending_longitude_moon(jce)
    ]
    let sum_psi = 0, sum_epsilon = 0;
    PE_TERMS.forEach((PE, i) => {
        let xy_term_sum = 0;

        for (let j = 0; j < data.length; j++) {
            xy_term_sum += data[j] * Y_TERMS[i][j];
        }

        sum_psi += (PE[0] + jce * PE[1]) * sin(rad * xy_term_sum);
        sum_epsilon += (PE[2] + jce * PE[3]) * cos(rad * xy_term_sum);
    })

    return [sum_psi / 36000000.0, sum_epsilon / 36000000.0];
}

function ecliptic_mean_obliquity(jme) {
    let u = jme / 10.0;

    return 84381.448 + u * (-4680.93 + u * (-1.55 + u * (1999.25 + u * (-51.38 + u * (-249.67 +
        u * (-39.05 + u * (7.12 + u * (27.87 + u * (5.79 + u * 2.45)))))))));
}
const ecliptic_true_obliquity = (delta_epsilon, epsilon0) => {
    return delta_epsilon + epsilon0 / 3600.0;
}

let L = earth_heliocentric_longitude(JME(2020, 5, 27, 12, 43))
let B = earth_heliocentric_latitude(JME(2020, 5, 27, 12, 43))
let auR = earth_radius_vector(JME(2020, 5, 27, 12, 43))
let Theta = geocentric_longitude(L)
let Beta = geocentric_latitude(B);


// mean sideral time in degrees
// https://de.wikipedia.org/wiki/Sternzeit#Sternzeit_in_Greenwich
export function mean_sideral_time(y, m, d, hh, mm, ss) {
    const t = T(y, m, d, 0, 0, 0),
        mjd = MJD(y, m, d, hh, mm, ss);
    return ((24110.54841 + t * (8640184.812866 + t * (0.093104 - t * 6.2e-6))
        + (mjd % 1) * 86400 * 1.00273790934) % 86400)
        / 3600 * 15;
}

// mean obliquity of the ecliptic in radians
// https://www.nrel.gov/docs/fy08osti/34302.pdf
export function epsilon(y, m, d, hh = 0, mm = 0, ss = 0) {
    const t = T(y, m, d, hh, mm, ss) / 100;
    return rad * (84381.448 + t * (-4680.93 + t * (-1.55 + t * (1999.25 + t * (-51.38 + t * (-249.67 +
        t * (-39.05 + t * (7.12 + t * (27.87 + t * (5.79 + t * 2.45)))))))))) / 3600;
}


const aberration_correction = (auR) => { return -20.4898 / (3600.0 * auR); }
const jce = JCE(2020, 5, 27, 12, 43);


const [latitude, longitude] = [51.3595977774021, -0.0988603855692327];
//log()
let [Delta_psi, Del_epsilon] = nutation_longitude_and_obliquity(jce)// -0.004976,-0.000121
let eps0 = epsilon(2020, 5, 27, 12, 43);
let eps = ecliptic_true_obliquity(Del_epsilon, eps0 * 3600);
//log(del_psi, del_epsilon, eps0*3600*deg /*84371.897721*/, eps*deg /*23.436517?? */) 
let Delta_tau = aberration_correction(auR)

const apparent_sun_longitude = (Theta, Delta_psi, Delta_tau) => { return Theta + Delta_psi + Delta_tau; }

/* function greenwich_mean_sidereal_time(y, m, d, hh = 0, mm = 0, ss = 0) {
    const j = JD(y, m, d, hh, mm, ss),
        t = T(y, m, d, hh, mm, ss);
    return range(280.46061837 + 360.98564736629 * (j - 2451545.0) +
        t * t * (0.000387933 - t / 38710000.0));
} */

function greenwich_mean_sidereal_time(jd, jc) {
    return range(280.46061837 + 360.98564736629 * (jd - 2451545.0) +
        jc * jc * (0.000387933 - jc / 38710000.0));
}

function greenwich_sidereal_time(Nu0, Delta_psi, Epsilon) {
    return Nu0 + Delta_psi * cos(rad * (Epsilon));
}
let Nu0 = greenwich_mean_sidereal_time(2020, 5, 27, 12, 43);
let Nu = greenwich_sidereal_time(Nu0, Delta_psi, eps * deg)
// log(Nu0, Nu, eps) // 76.284241,76.279675
//log(delta_tau) // -0.005616595805931507

function geocentric_right_ascension(Lambda, Epsilon, Beta) {
    let lamda_rad = rad * (Lambda);
    let epsilon_rad = rad * (Epsilon);

    return range(deg * (atan2(sin(lamda_rad) * cos(epsilon_rad) -
        tan(rad * (Beta)) * sin(epsilon_rad), cos(lamda_rad))));
}

let Lambda = apparent_sun_longitude(Theta, Delta_psi, Delta_tau); // 66.68463997838464 . 66.684564
let Ra = geocentric_right_ascension(Lambda, eps * deg, Beta); // 64.83985380190462 (Alpha)

function geocentric_declination(Beta, Epsilon, Lambda) {
    let beta_rad = rad * (Beta);
    let epsilon_rad = rad * (Epsilon);

    return deg * (asin(sin(beta_rad) * cos(epsilon_rad) +
        cos(beta_rad) * sin(epsilon_rad) * sin(rad * (Lambda))));
}

let Delta = geocentric_declination(Beta, eps * deg, Lambda); // 21.417165881054775

function observer_hour_angle(Nu, longitude, Ra) {
    return range(Nu + longitude - Ra);
}

let H = observer_hour_angle(Nu, longitude, Ra); // 11.340960581116008

function sun_equatorial_horizontal_parallax(auR) {
    return 8.794 / (3600.0 * auR);
}

let Xi = sun_equatorial_horizontal_parallax(auR); //0.002410582022145735

function right_ascension_parallax_and_topocentric_dec(latitude, elevation, Xi, H, Delta) {
    let lat_rad = rad * (latitude);
    let xi_rad = rad * (Xi);
    let h_rad = rad * (H);
    let delta_rad = rad * (Delta);
    let u = atan(0.99664719 * tan(lat_rad));
    let y = 0.99664719 * sin(u) + elevation * sin(lat_rad) / 6378140.0;
    let x = cos(u) + elevation * cos(lat_rad) / 6378140.0;

    let delta_alpha_rad = atan2(- x * sin(xi_rad) * sin(h_rad),
        cos(delta_rad) - x * sin(xi_rad) * cos(h_rad));

    let delta_prime = deg * (atan2((sin(delta_rad) - y * sin(xi_rad)) * cos(delta_alpha_rad),
        cos(delta_rad) - x * sin(xi_rad) * cos(h_rad)));

    let delta_alpha = deg * (delta_alpha_rad);
    return [delta_alpha, delta_prime]
}

let [Delta_alpha, Delta_prime] = right_ascension_parallax_and_topocentric_dec(latitude, 0/* elevation */, Xi, H, Delta) // -0.0003186175380131155 21.41596118399436

function topocentric_right_ascension(Ra, Delta_alpha) {
    return Ra + Delta_alpha;
}

let Alpha_prime = topocentric_right_ascension(Ra, Delta_alpha); // 64.83953518436661

function topocentric_local_hour_angle(H, Delta_alpha) {
    return H - Delta_alpha;
}

let H_prime = topocentric_local_hour_angle(H, Delta_alpha) // 11.341279198654021

function topocentric_elevation_angle(latitude, Delta_prime, H_prime) {
    let lat_rad = rad * (latitude);
    let delta_prime_rad = rad * (Delta_prime);

    return deg * (asin(sin(lat_rad) * sin(delta_prime_rad) +
        cos(lat_rad) * cos(delta_prime_rad) * cos(rad * (H_prime))));
}

let E0 = topocentric_elevation_angle(latitude, Delta_prime, H_prime); // 58.778036865420596

function atmospheric_refraction_correction(pressure, temperature, atmos_refract, E0) {
    let Del_e = 0;

    if (E0 >= -1 * (SUN_RADIUS + atmos_refract))
        Del_e = (pressure / 1010.0) * (283.0 / (273.0 + temperature)) *
            1.02 / (60.0 * tan(rad * (E0 + 10.3 / (E0 + 5.11))));

    return Del_e;
}

const pressure = 1000, temperature = 20, atmos_refract = 0.5667;
let Delta_e = atmospheric_refraction_correction(pressure, temperature, atmos_refract, E0) // 0.009791798142171966

function topocentric_elevation_angle_corrected(E0, Delta_e) { return E0 + Delta_e; }

let E = topocentric_elevation_angle_corrected(E0, Delta_e); // 58.78782866356277

function topocentric_zenith_angle(E) { return 90.0 - E; }

let Zenith = topocentric_zenith_angle(E) // 31.21217133643723

function topocentric_azimuth_angle_astro(H_prime, latitude, Delta_prime) {
    let h_prime_rad = rad * (H_prime);
    let lat_rad = rad * (latitude);

    return range(deg * (atan2(sin(h_prime_rad),
        cos(h_prime_rad) * sin(lat_rad) - tan(rad * (Delta_prime)) * cos(lat_rad))));
}

let Azimuth_astro = topocentric_azimuth_angle_astro(H_prime, latitude, Delta_prime) // Top. azimuth angle (westward from S) : 20.6821745142251

function topocentric_azimuth_angle(Azimuth_astro) { return range(Azimuth_astro + 180.0); }

let Azimuth = topocentric_azimuth_angle(Azimuth_astro) // Top. azimuth angle (eastward from N) : 200.68217451422512

function surface_incidence_angle(Zenith, Azimuth_astro, Azm_rotation, Slope) {
    let zenith_rad = rad * (Zenith);
    let slope_rad = rad * (Slope);

    return deg * (acos(cos(zenith_rad) * cos(slope_rad) +
        sin(slope_rad) * sin(zenith_rad) * cos(rad * (Azimuth_astro - Azm_rotation))));
}
// Surface azimuth rotation (measured from south to projection of surface normal on horizontal plane, negative east) [degrees]
// Surface slope (measured from the horizontal plane) [degrees]
const Azm_rotation = 180, Slope = 0;
let Incidence = surface_incidence_angle(Zenith, Azimuth_astro, Azm_rotation, Slope); // 31.21217133643723


function sun_mean_longitude(jme) {
    return range(280.4664567 + jme * (360007.6982779 + jme * (0.03032028 +
        jme * (1 / 49931.0 + jme * (-1 / 15300.0 + jme * (-1 / 2000000.0))))));
}
const jme = JME(2020, 5, 27, 12, 42);
let M = sun_mean_longitude(jme); // NO RESULTS AVAILABLE

function limit_minutes(minutes) {
    let limited = minutes;

    if (limited < -20.0) limited += 1440.0;
    else if (limited > 20.0) limited -= 1440.0;

    return limited;
}

// Equation of time
function eot(M, Alpha, Delta_psi, Epsilon) {
    return limit_minutes(4.0 * (M - 0.0057183 - Alpha + Delta_psi * cos(rad * (Epsilon))));
}

let Eot = eot(M, Ra, Delta_psi, eps * deg) // 2.7600946438849654

function approx_sun_transit_time(Ra0, longitude, Nu) {
    return (Ra0 - longitude - Nu) / 360.0;
}

let m0 = approx_sun_transit_time(Ra /*Ra0*/, longitude, Nu) // ???

const h0_prime = -1 * (SUN_RADIUS + atmos_refract);

const JD_MINUS = 0, JD_ZERO = 1, JD_PLUS = 2, JD_COUNT = 3;
const SUN_TRANSIT = 0, SUN_RISE = 1, SUN_SET = 2, SUN_COUNT = 3;

function sun_hour_angle_at_rise_set(latitude, delta_zero, h0_prime) {
    let H0 = -99999;
    let latitude_rad = latitude * rad;
    let delta_zero_rad = delta_zero * rad;

    let argument = (sin(rad * h0_prime) - sin(latitude_rad) * sin(delta_zero_rad))
        / (cos(latitude_rad) * cos(delta_zero_rad));
    if (Math.abs(argument) <= 1) H0 = range(deg * acos(argument), 180);
    return H0;
}
let H0 = sun_hour_angle_at_rise_set(latitude, Delta /*Delta0*/, h0_prime) // 121.04174815644211

let m_rts = [m0]

function approx_sun_rise_and_set(m_rts, H0) {
    const h0_dfrac = H0 / 360.0;

    m_rts[SUN_RISE] = range(m_rts[SUN_TRANSIT] - h0_dfrac, 1);
    m_rts[SUN_SET] = range(m_rts[SUN_TRANSIT] + h0_dfrac, 1);
    m_rts[SUN_TRANSIT] = range(m_rts[SUN_TRANSIT], 1);

    return m_rts;
}
m_rts = approx_sun_rise_and_set(m_rts, H0) // ??? 3.896418,11.960457,20.036836


log(m_rts)

// sun mean anomaly in degrees
// https://www.nrel.gov/docs/fy08osti/34302.pdf
export function sun_mean_anomaly(y, m, d, hh = 0, mm = 0, ss = 0) {
    const jce = JCE(y, m, d, hh, mm, ss);
    return (357.52772
        + 35999.050340 * jce
        - 0.0001603 * pow(jce, 2)
        - (pow(jce, 3) / 300000)
    ) % 360;
}

export function sun_right_ascension(l, b, e) { return atan2(sin(l) * cos(e) - tan(b) * sin(e), cos(l)); }
export function sun_declination(l, b, e) { return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l)); }

export function jdToTime(j) {
    const hh = ((j - 0.5) % 1) / (1 / 24),
        mm = (hh % 1) / (1 / 60),
        ss = (mm % 1) / (1 / 60);
    return [floor(hh), floor(mm), floor(ss)]
}

export function jdToDate(jd) {
    var l, n, i, j, k;

    l = jd + 68569;
    n = Math.floor(Math.floor(4 * l) / 146097);
    l = l - Math.floor((146097 * n + 3) / 4);
    i = Math.floor(4000 * (l + 1) / 1461001);
    l = l - Math.floor(1461 * i / 4) + 31;
    j = Math.floor(80 * l / 2447);
    k = l - Math.floor(2447 * j / 80);
    l = Math.floor(j / 11);
    j = j + 2 - 12 * l;
    i = 100 * (n - 49) + i + l;

    return [i, j, Math.round(k)]
}

export function jdToDateTime(jd) {
    return jdToDate(jd).concat(jdToTime(jd))
}
function ecliptical_longitude(M, C) { //  perihelion of the Earth
    return (M + 102.9373 + C + 180) % 360;
}
function eclipticLongitude(M) {

    var C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M)), // equation of center
        P = rad * 102.9372; // perihelion of the Earth

    return M + C + P + PI;
}

function equation_of_center(M) {
    return 1.9148 * sin(M)
        + 0.0200 * sin(2 * M)
        + 0.0003 * sin(3 * M);
}

function solar_transit2(lon, year, month, day, hour, mm, ss) {
    const M = sun_mean_anomaly(year, month, day, hour),
        C = equation_of_center(M),
        lambda = ecliptical_longitude(M, C),
        j = JD(year, month, day, hour, mm, ss),
        nx = j - J2000 - 0.0009 - ((lon * -1) / 360),
        jx = j + (trunc(nx) - nx);
    return jx + 0.0053 * sin(rad * M) - 0.0068 * sin(rad * (2 * lambda)); // (35)
}

function solar_transit(jd, M, lambda, lon) {
    let nx = jd - J2000 - 0.0009 - ((lon * -1) / 360); // (33)
    let jx = jd + (Math.trunc(nx) - nx) // (34)
    return jx + 0.0053 * sin(rad * M) - 0.0068 * sin(rad * (2 * lambda))// (35)
}

// Julian day,      Julian century,     Julian ephemeris day,   Julian ephemeris century,   Julian ephemeris millennium
// 2458997.029861,  0.204025,           2458997.030611,         0.204025,                   0.020403
/* let [y, m, d, hh, mm] = [2020, 5, 27, 12, 43];
log(JD(y, m, d, hh, mm), T(y, m, d, hh, mm), JDE(y, m, d, hh, mm), JCE(y, m, d, hh, mm), JME(y, m, d, hh, mm)) */

const [year, month, day, hh, mm, ss] = [2020, 5, 27, 12, 43, 0];

function julian_century(jd) {
    return (jd - 2451545.0) / 36525.0;
}

function julian_ephemeris_day(jd, delta_t) {
    return jd + delta_t / 86400.0;
}

function julian_ephemeris_century(jde) {
    return (jde - 2451545.0) / 36525.0;
}

function julian_ephemeris_millennium(jce) {
    return (jce / 10.0);
}
function calculate_geocentric_sun_right_ascension_and_declination(spa) {

   // spa.jd = JD(year, month, day, hh, mm, ss);
    spa.jc = julian_century(spa.jd);//JC(year, month, day, hh, mm, ss);
    spa.jde = julian_ephemeris_day(spa.jd, spa.delta_t)//JDE(year, month, day, hh, mm, ss);
    spa.jce = julian_ephemeris_century(spa.jde)//JCE(year, month, day, hh, mm, ss);
    spa.jme = julian_ephemeris_millennium(spa.jce)//JME(year, month, day, hh, mm, ss);
    spa.l = earth_heliocentric_longitude(spa.jme);
    spa.b = earth_heliocentric_latitude(spa.jme);
    spa.r = earth_radius_vector(spa.jme);
    spa.theta = geocentric_longitude(spa.l);
    spa.beta = geocentric_latitude(spa.b);
    [spa.del_psi, spa.del_epsilon] = nutation_longitude_and_obliquity(spa.jce)

    spa.epsilon0 = ecliptic_mean_obliquity(spa.jme);
    spa.epsilon = ecliptic_true_obliquity(spa.del_epsilon, spa.epsilon0);

    spa.del_tau = aberration_correction(spa.r);
    spa.lamda = apparent_sun_longitude(spa.theta, spa.del_psi, spa.del_tau);
    //spa.nu0 = greenwich_mean_sidereal_time(year, month, day, hh, mm, ss)//spa.jd, spa.jc); // TODO
    spa.nu0 = greenwich_mean_sidereal_time(spa.jd, spa.jc);
    spa.nu = greenwich_sidereal_time(spa.nu0, spa.del_psi, spa.epsilon);

    spa.alpha = geocentric_right_ascension(spa.lamda, spa.epsilon, spa.beta);
    spa.delta = geocentric_declination(spa.beta, spa.epsilon, spa.lamda);
}

function spa_calculate(spa) {
    spa.jd = JD(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.second, spa.delta_ut1, spa.timezone);

    calculate_geocentric_sun_right_ascension_and_declination(spa);

    spa.h = observer_hour_angle(spa.nu, spa.longitude, spa.alpha);
    spa.xi = sun_equatorial_horizontal_parallax(spa.r);

    [spa.del_alpha, spa.delta_prime] = right_ascension_parallax_and_topocentric_dec(spa.latitude, spa.elevation, spa.xi, spa.h, spa.delta);

    spa.alpha_prime = topocentric_right_ascension(spa.alpha, spa.del_alpha);
    spa.h_prime = topocentric_local_hour_angle(spa.h, spa.del_alpha);

    spa.e0 = topocentric_elevation_angle(spa.latitude, spa.delta_prime, spa.h_prime);
    spa.del_e = atmospheric_refraction_correction(spa.pressure, spa.temperature, spa.atmos_refract, spa.e0);
    spa.e = topocentric_elevation_angle_corrected(spa.e0, spa.del_e);

    spa.zenith = topocentric_zenith_angle(spa.e);
    spa.azimuth_astro = topocentric_azimuth_angle_astro(spa.h_prime, spa.latitude, spa.delta_prime);
    spa.azimuth = topocentric_azimuth_angle(spa.azimuth_astro);

    calculate_eot_and_sun_rise_transit_set(spa);
}


function rts_alpha_delta_prime(ad, n) {
    let a = ad[JD_ZERO] - ad[JD_MINUS];
    let b = ad[JD_PLUS] - ad[JD_ZERO];

    if (abs(a) >= 2.0) a = range(a, 1);
    if (abs(b) >= 2.0) b = range(b, 1);

    return ad[JD_ZERO] + n * (a + b + (b - a) * n) / 2.0;
}

function rts_sun_altitude(latitude, delta_prime, h_prime) {
    let latitude_rad = rad * (latitude);
    let delta_prime_rad = rad * (delta_prime);

    return deg * (asin(sin(latitude_rad) * sin(delta_prime_rad) +
        cos(latitude_rad) * cos(delta_prime_rad) * cos(rad * (h_prime))));
}

function dayfrac_to_local_hr(dayfrac, timezone) {
    return 24.0 * range(dayfrac + timezone / 24.0, 1);
}

function sun_rise_and_set(m_rts, h_rts, delta_prime, latitude, h_prime, h0_prime, sun) {
    return m_rts[sun] + (h_rts[sun] - h0_prime) /
        (360.0 * cos(rad * (delta_prime[sun])) * cos(rad * (latitude)) * sin(rad * (h_prime[sun])));
}

function calculate_eot_and_sun_rise_transit_set(spa) {
    // spa_data sun_rts;
    // double nu, m, h0, n;
    let alpha = [], delta = [],
        m_rts = [], nu_rts = [], h_rts = [],                  // double alpha[JD_COUNT], delta[JD_COUNT];
        // double m_rts[SUN_COUNT], nu_rts[SUN_COUNT], h_rts[SUN_COUNT];
        // double alpha_prime[SUN_COUNT], delta_prime[SUN_COUNT], h_prime[SUN_COUNT];
        alpha_prime = [], delta_prime = [], h_prime = [];
    let h0_prime = -1 * (SUN_RADIUS + spa.atmos_refract);   // double h0_prime = -1*(SUN_RADIUS + spa.atmos_refract);
    // int i;

    let sun_rts = { ...spa },
        m = sun_mean_longitude(spa.jme);
    spa.eot = eot(m, spa.alpha, spa.del_psi, spa.epsilon);

    sun_rts.hour = sun_rts.minute = sun_rts.second = 0;
    sun_rts.delta_ut1 = sun_rts.timezone = 0.0;

    sun_rts.jd = JD(sun_rts.year, sun_rts.month, sun_rts.day, sun_rts.hour,
        sun_rts.minute, sun_rts.second, sun_rts.delta_ut1, sun_rts.timezone);

    calculate_geocentric_sun_right_ascension_and_declination(sun_rts);

    let nu = sun_rts.nu;

    sun_rts.delta_t = 0;
    sun_rts.jd--;

    for (let i = 0; i < JD_COUNT; i++) {
        calculate_geocentric_sun_right_ascension_and_declination(sun_rts);
        alpha[i] = sun_rts.alpha;
        delta[i] = sun_rts.delta;
        sun_rts.jd++;
    }
    
    m_rts[SUN_TRANSIT] = approx_sun_transit_time(alpha[JD_ZERO], spa.longitude, nu);
    let h0 = sun_hour_angle_at_rise_set(spa.latitude, delta[JD_ZERO], h0_prime);

    if (h0 >= 0) {
        approx_sun_rise_and_set(m_rts, h0);

        for (let i = 0; i < SUN_COUNT; i++) {
            nu_rts[i] = nu + 360.985647 * m_rts[i];

            let n = m_rts[i] + spa.delta_t / 86400.0;
            alpha_prime[i] = rts_alpha_delta_prime(alpha, n);
            delta_prime[i] = rts_alpha_delta_prime(delta, n);

            h_prime[i] = range(nu_rts[i] + spa.longitude - alpha_prime[i], 180, -180);

            h_rts[i] = rts_sun_altitude(spa.latitude, delta_prime[i], h_prime[i]);
        }
log(h_prime,h_rts, nu_rts,m_rts )
        spa.srha = h_prime[SUN_RISE];
        spa.ssha = h_prime[SUN_SET];
        spa.sta = h_rts[SUN_TRANSIT];

        spa.suntransit = dayfrac_to_local_hr(m_rts[SUN_TRANSIT] - h_prime[SUN_TRANSIT] / 360.0, spa.timezone);

        spa.sunrise = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,
            spa.latitude, h_prime, h0_prime, SUN_RISE), spa.timezone);

        spa.sunset = dayfrac_to_local_hr(sun_rise_and_set(m_rts, h_rts, delta_prime,
            spa.latitude, h_prime, h0_prime, SUN_SET), spa.timezone);
    }
}

function main() {
    let spa = {
        year: 2020,
        month: 5,
        day: 27,
        hour: 12,
        minute: 43,
        second: 0,
        delta_ut1: 0,
        delta_t: 64.8,
        timezone: 0,
        longitude: longitude,
        latitude: latitude,
        elevation: 0,
        pressure: 1000,
        temperature: 20,
        slope: 0,
        azm_rotation: 0,
        atmos_refract: 0.5667
    };

    spa_calculate(spa);


/*         for (var i in spa)
            log(i, spa[i]) */
    log(spa.sunrise,spa.suntransit,spa.sunset) // 3.896418,11.960457,20.036836
    log(hhmmss(spa.sunrise).join(':'),hhmmss(spa.suntransit).join(':'),hhmmss(spa.sunset).join(':'))
}
main()