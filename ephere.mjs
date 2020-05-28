/* EPHERE.MJS */

const PI = Math.PI,
    J2000 = 2451545.0,
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
const log = console.log;

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
const JDE = (y, m, d, hh, mm, ss) => { return JD(y, m, d, hh, mm, ss) + (deltaT(y) / 86400); }
const JCE = (y, m, d, hh, mm, ss) => { return (JDE(y, m, d, hh, mm, ss) - J2000) / 36525; }

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
    return rad *(84381.448
        - 4680.93 * t
        - 1.55 * pow(t, 2)
        + 1999.25 * pow(t, 3)
        - 51.38 * pow(t, 4)
        - 249.67 * pow(t, 5)
        - 39.05 * pow(t, 6)
        + 7.12 * pow(t, 7)
        + 27.87 * pow(t, 8)
        + 5.79 * pow(t, 9)
        + 2.45 * pow(t, 10)) / 3600;
}

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

