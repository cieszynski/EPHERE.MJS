/* sun.mjs */
import { CelestialObject } from './ephere.mjs'

export default class Sun extends CelestialObject {
    xyz = (t, origin) => {
        return [0, 0, 0].map((n, i) => n - origin[i]);
    }
}