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
//const range = (x, max = 360) => { return max * ((x / max) % 1) + (x < 0) * max; }
const range = (x, max = 360, min = 0) => {
    return ((max - min) * frac(x / (max - min)))
        + ((max - min) * frac(x / (max - min)) < min) * (max - min)
        - ((max - min) * frac(x / (max - min)) > max) * (max - min);
}
