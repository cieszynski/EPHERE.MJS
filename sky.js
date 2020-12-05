/* sky.js */
"use strict"

class Sun extends CelestialObject {
    constructor() {
        super();
        this.h0 = -0.833;
    }

    xyz(t, origin) {
        const v = [0, 0, 0].map(function (n, i) { return n - origin[i] });
        return {
            x: v[0],
            y: v[1],
            z: v[2]
        }
    }
}

class Projection {

}

class Lambert extends Projection {

    azel2xy(az, al, w, h) {
        const cosaz = Math.cos((az - Math.PI));
        const sinaz = Math.sin((az - Math.PI));
        const sinel = Math.sin(al);
        const cosel = Math.cos(al);
        const k = Math.sqrt(2 / (1 + cosel * cosaz));
        return { x: (w / 2 + 0.6 * h * k * cosel * sinaz), y: (h - 0.6 * h * k * (sinel)), al: al };
    }
}

class Polar extends Projection {

    azel2xy(az, al, w, h) {
        var radius = h / 2;
        var r = radius * ((Math.PI / 2) - al) / (Math.PI / 2);
        return { x: (w / 2 - r * Math.sin(az)), y: (radius - r * Math.cos(az)), al: al };
    }
}

class EquiRectangular extends Projection {

    azel2xy(az, al, w, h) {
        while (az < 0) az += 2 * Math.PI;
        az = (az) % (Math.PI * 2);
        return { x: (((az - Math.PI) / (Math.PI / 2)) * h + w / 2), y: (h - (al / (Math.PI / 2)) * h), al: al };
    }
}

class StereoGraphic extends Projection {

    azel2xy(az, al, w, h) {
        var f = 0.5
        var sinel1 = 0;
        var cosel1 = 1;
        var cosaz = Math.cos((az - Math.PI));
        var sinaz = Math.sin((az - Math.PI));
        var sinel = Math.sin(al);
        var cosel = Math.cos(al);
        var k = 2 / (1 + sinel1 * sinel + cosel1 * cosel * cosaz);
        return { x: (w / 2 + f * k * h * cosel * sinaz), y: (h - f * k * h * (cosel1 * sinel - sinel1 * cosel * cosaz)), al: al };
    }
}

class Sky {

    constructor(canvas, opts) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d', { alpha: false });
        this.width = opts.width || 640;
        this.height = opts.height || 320;

        this.date = opts.date || new Date();
        this.jd = JD(
            this.date.getUTCFullYear(),
            this.date.getUTCMonth() + 1,
            this.date.getUTCDate(),
            this.date.getUTCHours(),
            this.date.getUTCMinutes(),
            this.date.getUTCSeconds());
        this.lat = opts.lat || 51.477778;   // Greenwich Observatory
        this.lon = opts.lon || -0.001389;   // Greenwich Observatory

        this.projection = opts.projection || new Lambert();
        this.cardinalpoints = ['N', 'E', 'S', 'W'];
        this.cardinalpoints_font = opts.cardinalpoints_font || 'sans-serif';
        this.cardinalpoints_fontsize = opts.cardinalpoints_fontsize || 10;
        this.cardinalpoints_color = opts.cardinalpoints_color || 'gray';
        this.cardinalpoints_bottom = opts.cardinalpoints_bottom || 15;

        this.constellations = opts.constellations || [];
        this.constellations_linewidth = opts.constellations_linewidth || 0.75;
        this.constellations_linecolor = opts.constellations_linecolor || 'gray';
        this.constellations_fontcolor = opts.constellations_fontcolor || 'lightgray';
        this.constellations_fontsize = opts.constellations_fontsize || 12;
        this.constellations_font = opts.constellations_font || 'sans-serif';

        this.planets = opts.planets || [];
        this.planets_size = opts.planets_size || 5;
        this.planets_color = opts.planets_color || 'white';
        this.planets_font = opts.planets_font || 'sans-serif';
        this.planets_fontsize = opts.planets_fontsize || 14;
        this.stars = opts.stars || [];
        this.scalestars = opts.scalestars || 1;
        this.i18n = opts.i18n || [];

        this.x = 0;
        this.y = 0;
        this.az_off = 0;
        this.dragging = false;

        this.canvas.ontouchstart = this.ontouchstart.bind(this);
        this.canvas.ontouchend = this.ontouchend.bind(this);
        this.canvas.ontouchmove = this.ontouchmove.bind(this);
        this.sun = new Sun();

    }

    set width(value) { this.canvas.width = value; }
    get width() { return this.canvas.width; }

    set height(value) { this.canvas.height = value; }
    get height() { return this.canvas.height; }

    ontouchstart(e) {
        this.x = e.touches[0].pageX;
        this.y = e.touches[0].pageY;
        this.dragging = true;
    }

    ontouchend(e) { this.dragging = false; }

    ontouchmove(e) {
        if (this.move)
            return
        const x = e.touches[0].pageX;
        const y = e.touches[0].pageY;

        if (this.dragging) {
            this.az_off += (this.x - x) / 4;
            this.az_off = this.az_off % 360;
            this.x = x;
            this.y = y;
            this.canvas.ontouchmove = null;
            requestAnimationFrame(this.draw.bind(this))
        }
    }

    font(fs, f) {
        const m = Math.min(this.width, this.height);
        //fs *= (m < 800) ? ((m < 700) ? ((m < 600) ? ((m < 500) ? (m < 400) ? 1 : 1.2 : 1.4) : 1.6) : 1.8) : 2;
        return `${fs}px ${f}`;
    }

    skygradient() {

        const gradient = this.ctx.createLinearGradient(0, 0, 0, this.height);
        gradient.addColorStop(0.0, 'rgba(0,30,50,0.1)');
        gradient.addColorStop(0.7, 'rgba(0,30,50,0.35)');
        gradient.addColorStop(1, 'rgba(0,50,80,0.6)');

        return gradient;
    }

    draw() {
        this.ctx.clearRect(0, 0, this.width, this.height);
        this.ctx.fillStyle = 'black';
        this.ctx.fillRect(0, 0, this.width, this.height);
        this.ctx.fill();

        this.ctx.beginPath();
        this.ctx.fillStyle = this.skygradient();
        this.ctx.fillRect(0, 0, this.width, this.height);
        this.ctx.closePath();

        this.drawCardinalPoints()
        this.drawStars();
        this.drawPlanets();
        this.drawConstellations();
        this.canvas.ontouchmove = this.ontouchmove.bind(this);
    }

    drawCardinalPoints() {
        const cp = this.cardinalpoints;
        const cpa = cp.map(function (value, idx, arr) {
            return 360 / arr.length * idx;
        });
        const cpf = this.cardinalpoints_font;
        const cpfs = this.cardinalpoints_fontsize;
        const cpc = this.cardinalpoints_color;
        const cpb = this.cardinalpoints_bottom;

        this.ctx.beginPath();
        this.ctx.fillStyle = cpc;
        this.ctx.font = this.font(cpfs, cpf);

        for (let i = 0; i < cp.length; i++) {
            const ang = (cpa[i] - this.az_off) * rad;
            const m = this.ctx.measureText(cp[i]);
            const r = (m.width > cpfs) ? m.width / 2 : cpfs / 2;
            const pos = this.projection.azel2xy(ang, 0, this.width, this.height);

            let x = isFinite(pos.x) ? pos.x - r : 0;
            let y = isFinite(pos.y) ? pos.y - cpb / 2 : 0;

            if (x < 0 || x > this.width - r) x = -r;
            if (x > 0) this.ctx.fillText(this.i18n[cp[i]], x, y);
        }

        this.ctx.fill();
    }

    drawStars() {
        this.ctx.beginPath();
        this.ctx.fillStyle = 'lightgray'

        for (let i = 0; i < this.stars.length; i++) {
            const pos = this.radec2xy(this.stars[i][2], this.stars[i][3]);

            if (!this.isvisible(pos))
                continue;

            const mag = this.stars[i][1];
            const d = this.scalestars * Math.max(3 - mag / 2.1, 0.5);
            this.ctx.moveTo(pos.x + d, pos.y);
            this.ctx.arc(pos.x, pos.y, d, 0, Math.PI * 2, true);
        }
        this.ctx.fill();
    }

    drawConstellations() {
        this.ctx.fillStyle = this.constellations_textcolor;
        this.ctx.strokeStyle = this.constellations_linecolor;
        this.ctx.lineWidth = this.constellations_linewidth;
        this.ctx.font = this.font(this.constellations_fontsize, this.constellations_font);
        this.ctx.beginPath();

        const maxl = this.height / 3;
        const hipparcos = {}

        for (let c = 0; c < this.constellations.length; c++) {
            for (let l = 3; l < this.constellations[c].length; l += 2) {
                let a = -1;
                let b = -1;
                let idx1 = '' + this.constellations[c][l] + '';
                let idx2 = '' + this.constellations[c][l + 1] + '';
                if (!hipparcos[idx1]) {
                    for (let s = 0; s < this.stars.length; s++) {
                        if (this.stars[s][0] == this.constellations[c][l]) {
                            hipparcos[idx1] = s;
                            break;
                        }
                    }
                }
                if (!hipparcos[idx2]) {
                    for (let s = 0; s < this.stars.length; s++) {
                        if (this.stars[s][0] == this.constellations[c][l + 1]) {
                            hipparcos[idx2] = s;
                            break;
                        }
                    }
                }
                a = hipparcos[idx1];
                b = hipparcos[idx2];
                if (a >= 0 && b >= 0 && a < this.stars.length && b < this.stars.length) {
                    let posa = this.radec2xy(this.stars[a][2], this.stars[a][3]);
                    let posb = this.radec2xy(this.stars[b][2], this.stars[b][3]);

                    if (!this.isvisible(posa) && !this.isvisible(posb))
                        continue;
                    //  if (!isPointBad(posa) && !isPointBad(posb)) {
                    // Basic error checking: i18n behind us often have very long lines so we'll zap them
                    if (Math.abs(posa.x - posb.x) < maxl && Math.abs(posa.y - posb.y) < maxl) {
                        this.ctx.moveTo(posa.x, posa.y);
                        this.ctx.lineTo(posb.x, posb.y);
                    }
                    //  }
                    //}
                }
            }

            let pos = this.radec2xy(this.constellations[c][1] * rad, this.constellations[c][2] * rad);
            let i18n = this.i18n[this.constellations[c][0]];
            let xoff = -this.ctx.measureText(i18n).width / 2;
            this.ctx.fillStyle = this.constellations_fontcolor;
            this.ctx.fillText(i18n, pos.x + xoff, pos.y - this.constellations_fontsize / 2);

        }
        this.ctx.fill();
        this.ctx.stroke();
    }

    drawPlanets() {
        for (let p = 0; p < this.planets.length; p++) {
            let d = this.planets_size;
            const name = this.planets[p].constructor.name;
            const label = this.i18n[name];
            const rradec = this.planets[p].rradec(this.jd);
            const pos = this.radec2xy(rradec.ra, rradec.dec);

            const x = pos.x;
            const y = pos.y;

            this.ctx.font = this.font(this.planets_fontsize, this.planets_font);
            this.ctx.beginPath();

            switch (name) {
                case "Moon":
                    d = 20;
                    this.ctx.fillStyle = this.planets_color;
                    this.ctx.strokeStyle = this.planets_color;
                    this.ctx.moveTo(x + d, y + d);
                    this.ctx.arc(x, y, d, 0, Math.PI * 2, true);
                    break;
                case "Sun":
                    //console.log(pos.al)
                    let D1 = (pos.al > 0) ? 30 : 0;
                    let D2 = (pos.al > 0) ? 100 : 0;
                    //this.ctx.fillStyle = '#f8f9f3';
                    this.ctx.moveTo(x + d, y + d);
                    this.ctx.arc(x, y, d, 0, Math.PI * 2, true);
                    const lineargradient = this.ctx.createRadialGradient(x, y, D1, x, y, D2)
                    lineargradient.addColorStop(0, '#f8f9f3');
                    lineargradient.addColorStop(1, 'rgba(255,255,255,0)');
                    this.ctx.fillStyle = lineargradient
                    this.ctx.moveTo(x + d, y + d);
                    this.ctx.arc(x, y, d, 0, Math.PI * 2, true);
                    break;
                default:
                    this.ctx.fillStyle = this.planets_color;
                    this.ctx.moveTo(x + d, y + d);
                    this.ctx.arc(x, y, d, 0, Math.PI * 2, true);
                    this.ctx.fillText(label, x + d, y - (d + 2));
            }
            this.ctx.fill();
        }
    }

    radec2xy(ra, dec) {
        const coords = this.sun.azal(this.jd, this.lat * rad, this.lon * rad, ra, dec)
        var pos = this.projection.azel2xy(coords.az - (this.az_off * rad), coords.al, this.width, this.height);
        return { x: Math.round(pos.x), y: Math.round(pos.y), az: coords.az * deg, al: coords.al * deg };
    };

    isvisible(pos) {
        return (pos !== undefined && pos.x >= 0 && pos.x <= this.width && pos.al > 0);
    }
}