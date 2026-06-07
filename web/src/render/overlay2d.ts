/* Canvas2D overlay drawn on top of the WebGL field canvas: particle trails,
 * markers, velocity arrows, and selection highlight.  Sim coordinates have
 * y=0 at the bottom (matching the WebGL particle pass); Canvas2D is y-down,
 * so we flip y in the mapping. */

import type { Particle } from '../sim/world';
import type { VectorData } from '../sim/fieldViews';
import type { TickMode, ForceArrows } from '../sim/scenario';

export const VEL_PX_PER_C = 420;   // arrow pixels at |v| = c
const TRAIL_MAX = 1200;      // points kept per particle (long enough to show a precessing rosette)
const TRAIL_STRIDE = 4;      // [x, y, t, tau] per trail point
const CLOCK_PERIOD = 120.0;  // time units per full hand revolution (gentle sweep)
const FORCE_REF_PX = 70;     // px length of the strongest displayed force arrow

/* Force-component arrow colors (§17): grav orange, EM blue, spin green,
 * total white. */
const FORCE_COLORS = {
    grav: '#ff9030', em: '#4aa3ff', spin: '#46d66a', total: '#f0f0f0',
} as const;

/* Per-particle visualization flags passed to render (parallel to particles). */
export interface ParticleDisplay {
    ticks: TickMode;
    clock: boolean;
    rate: number;   // instantaneous dτ/dt (for the clock readout)
}

export interface Trails {
    /* trails[i] is a flat ring of [x0,y0,t0,tau0, x1,y1,t1,tau1, ...]. */
    pts: number[][];
    push(particles: Particle[], t: number): void;
    reset(n: number): void;
}

export function createTrails(): Trails {
    const t: Trails = {
        pts: [],
        reset(n: number) { t.pts = Array.from({ length: n }, () => []); },
        push(particles: Particle[], time: number) {
            if (t.pts.length !== particles.length) t.reset(particles.length);
            for (let i = 0; i < particles.length; i++) {
                const arr = t.pts[i];
                arr.push(particles[i].x, particles[i].y, time, particles[i].properTime);
                const cap = TRAIL_MAX * TRAIL_STRIDE;
                if (arr.length > cap) arr.splice(0, arr.length - cap);
            }
        },
    };
    return t;
}

/* "Nice" number (1, 2, 5 × 10^k) at or below x -- for auto tick spacing. */
function niceInterval(x: number): number {
    if (!(x > 0)) return 1;
    const p = Math.pow(10, Math.floor(Math.log10(x)));
    const f = x / p;
    return (f >= 5 ? 5 : f >= 2 ? 2 : 1) * p;
}

export class Overlay2D {
    private ctx: CanvasRenderingContext2D;
    private cw: number;
    private ch: number;
    readonly W: number;
    readonly H: number;

    constructor(canvas: HTMLCanvasElement, simW: number, simH: number) {
        const ctx = canvas.getContext('2d');
        if (!ctx) throw new Error('2D context unavailable');
        this.ctx = ctx;
        this.cw = canvas.width;
        this.ch = canvas.height;
        this.W = simW; this.H = simH;
    }

    private sx(x: number): number { return (x / this.W) * this.cw; }
    private sy(y: number): number { return this.ch - (y / this.H) * this.ch; }

    /* nearest particle to a canvas pixel, or -1 if none within `tol` px. */
    pick(px: number, py: number, particles: Particle[], tol = 14): number {
        let best = -1, bestD = tol * tol;
        for (const p of particles) {
            const dx = this.sx(p.x) - px, dy = this.sy(p.y) - py;
            const d = dx * dx + dy * dy;
            if (d < bestD) { bestD = d; best = p.index; }
        }
        return best;
    }

    render(particles: Particle[], trails: Trails | null,
           opts: { showTrails: boolean; showVelocity: boolean;
                   selectedParticle: number; selectedBody?: number;
                   vectors?: { data: VectorData; spacing: number } | null;
                   display?: ParticleDisplay[]; tickInterval?: number; time?: number;
                   forces?: ForceArrows;
                   bodies?: { x: number; y: number; r: number }[] }): void {
        const ctx = this.ctx;
        const S = TRAIL_STRIDE;
        ctx.clearRect(0, 0, this.cw, this.ch);

        /* Vector-field arrows (under the trails/particles). */
        if (opts.vectors) this.drawVectors(opts.vectors.data, opts.vectors.spacing);

        /* Trails: fade with age (oldest transparent). */
        if (opts.showTrails && trails) {
            ctx.lineWidth = 1.5;
            for (let i = 0; i < trails.pts.length; i++) {
                const a = trails.pts[i];
                const segs = a.length / S;
                if (segs < 2) continue;
                for (let k = 1; k < segs; k++) {
                    const alpha = k / segs;
                    ctx.strokeStyle = `rgba(120,180,255,${(0.05 + 0.5 * alpha).toFixed(3)})`;
                    ctx.beginPath();
                    ctx.moveTo(this.sx(a[S * (k - 1)]), this.sy(a[S * (k - 1) + 1]));
                    ctx.lineTo(this.sx(a[S * k]), this.sy(a[S * k + 1]));
                    ctx.stroke();
                }
            }
        }

        /* Background compact body: a hollow bright-red circle at its position,
         * sized to its softening length (the body's effective core radius). */
        if (opts.bodies) {
            ctx.lineWidth = 2;
            for (let bi = 0; bi < opts.bodies.length; bi++) {
                const b = opts.bodies[bi];
                const rPx = Math.max(6, b.r * (this.cw / this.W));
                const bx = this.sx(b.x), by = this.sy(b.y);
                ctx.strokeStyle = '#ff2020';
                ctx.beginPath();
                ctx.arc(bx, by, rPx, 0, Math.PI * 2);
                ctx.stroke();
                /* Selected background body: amber highlight ring + center dot. */
                if (bi === opts.selectedBody) {
                    ctx.strokeStyle = '#ffe070';
                    ctx.beginPath(); ctx.arc(bx, by, rPx + 4, 0, Math.PI * 2); ctx.stroke();
                    ctx.fillStyle = '#ffe070';
                    ctx.beginPath(); ctx.arc(bx, by, 3, 0, Math.PI * 2); ctx.fill();
                }
            }
        }

        /* Coordinate/proper-time ticks along each selected trail. */
        if (trails && opts.display) {
            for (let i = 0; i < trails.pts.length; i++) {
                const d = opts.display[i];
                if (!d || d.ticks === 'none') continue;
                this.drawTrackTicks(trails.pts[i], d.ticks, opts.tickInterval ?? 0);
            }
        }

        for (const p of particles) {
            const x = this.sx(p.x), y = this.sy(p.y);
            const disp = opts.display?.[p.index];

            /* Velocity arrow (v = p / (gamma m), c=1). */
            if (opts.showVelocity && p.mass > 0) {
                const pmag2 = p.px * p.px + p.py * p.py;
                const gamma = Math.sqrt(1 + pmag2 / (p.mass * p.mass));
                const vx = p.px / (gamma * p.mass);
                const vy = p.py / (gamma * p.mass);
                const ex = x + vx * VEL_PX_PER_C;
                /* vy is +up in sim, -down in canvas */
                const ey = y - vy * VEL_PX_PER_C;
                const speed = Math.hypot(vx, vy);
                /* Violet, distinct from the force-arrow palette; brighter when
                 * relativistic (|v| > 0.5 c). */
                ctx.strokeStyle = speed > 0.5 ? '#e0b0ff' : '#9b6cff';
                ctx.lineWidth = 2;
                ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(ex, ey); ctx.stroke();
                this.arrowHead(ex, ey, Math.atan2(ey - y, ex - x), ctx.strokeStyle);
            }

            /* Gyroscope spin indicator: a "hand" through the particle at angle
             * phi_spin.  As the spin precesses (Lense-Thirring), the hand
             * rotates -- frame-dragging made literally visible.  A dot marks
             * the +phi tip so the 180-deg ambiguity is broken. */
            if (p.spin !== 0) {
                const L = 20;
                const hx = Math.cos(p.phiSpin) * L;
                const hy = -Math.sin(p.phiSpin) * L;   // sim y is +up, canvas +down
                ctx.strokeStyle = '#66ffcc';
                ctx.lineWidth = 2;
                ctx.beginPath();
                ctx.moveTo(x - hx, y - hy); ctx.lineTo(x + hx, y + hy); ctx.stroke();
                ctx.fillStyle = '#66ffcc';
                ctx.beginPath(); ctx.arc(x + hx, y + hy, 3, 0, Math.PI * 2); ctx.fill();
            }

            /* Marker. */
            const sel = p.index === opts.selectedParticle;
            ctx.beginPath();
            ctx.arc(x, y, sel ? 7 : 5, 0, Math.PI * 2);
            ctx.fillStyle = sel ? '#ffe070' : '#ffffff';
            ctx.fill();
            if (sel) {
                ctx.strokeStyle = '#ffe070'; ctx.lineWidth = 1.5;
                ctx.beginPath(); ctx.arc(x, y, 12, 0, Math.PI * 2); ctx.stroke();
            }
            ctx.fillStyle = '#888'; ctx.font = '10px ui-monospace, monospace';
            ctx.fillText(`p${p.index}`, x + 9, y - 9);

            /* Proper-time clock: coordinate-time hand (gray) vs proper-time
             * hand (cyan); the angular lag between them IS the accumulated time
             * dilation t-τ.  Digital readout shows τ, t-τ, and dτ/dt. */
            if (disp?.clock) this.drawClock(x + 26, y - 26, opts.time ?? p.properTime, p.properTime, disp.rate);
        }

        /* Force-component arrows on top of the markers. */
        if (opts.forces) this.drawForces(particles, opts.forces);
    }

    /* Per-particle force arrows (§17): gravitational / EM / spin-gradient /
     * total.  All enabled components across all particles share one auto-scale
     * so their lengths are directly comparable; the strongest spans FORCE_REF_PX.
     * Components stack from the same origin (the particle center). */
    private drawForces(particles: Particle[], f: ForceArrows): void {
        const ctx = this.ctx;
        const comps: { on: boolean; gx: (p: Particle) => number; gy: (p: Particle) => number; c: string }[] = [
            { on: f.grav,  gx: (p) => p.fgravX, gy: (p) => p.fgravY, c: FORCE_COLORS.grav },
            { on: f.em,    gx: (p) => p.femX,   gy: (p) => p.femY,   c: FORCE_COLORS.em },
            { on: f.spin,  gx: (p) => p.fspinX, gy: (p) => p.fspinY, c: FORCE_COLORS.spin },
            { on: f.total, gx: (p) => p.ftotX,  gy: (p) => p.ftotY,  c: FORCE_COLORS.total },
        ];
        let max = 0;
        for (const p of particles)
            for (const cm of comps)
                if (cm.on) max = Math.max(max, Math.hypot(cm.gx(p), cm.gy(p)));
        if (max <= 0) return;
        const scale = FORCE_REF_PX / max;

        ctx.lineWidth = 2;
        for (const p of particles) {
            const x = this.sx(p.x), y = this.sy(p.y);
            for (const cm of comps) {
                if (!cm.on) continue;
                const fx = cm.gx(p), fy = cm.gy(p);
                if (fx === 0 && fy === 0) continue;
                const ex = x + fx * scale, ey = y - fy * scale;   // sim y is +up
                ctx.strokeStyle = cm.c;
                ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(ex, ey); ctx.stroke();
                this.arrowHead(ex, ey, Math.atan2(ey - y, ex - x), cm.c, 7);
            }
        }
    }

    /* Place tick markers along a trail at equal coordinate-time (circles) and/or
     * proper-time (triangles) intervals.  Ticks sit at absolute multiples of Δ
     * so they stay pinned to the trajectory as the trail scrolls.  When Δ<=0,
     * Δ is auto-chosen from the trail's coordinate-time span (~12 ticks). */
    private drawTrackTicks(a: number[], mode: TickMode, interval: number): void {
        const S = TRAIL_STRIDE;
        const n = a.length / S;
        if (n < 2) return;
        let dv = interval;
        if (!(dv > 0)) {
            const span = a[(n - 1) * S + 2] - a[2];   // coordinate-time span
            dv = niceInterval(span / 12);
        }
        if (!(dv > 0)) return;
        if (mode === 'coordinate' || mode === 'both') this.ticksFor(a, 2, dv, false);
        if (mode === 'proper' || mode === 'both') this.ticksFor(a, 3, dv, true);
    }

    /* Walk trail segments emitting a marker wherever the value at stride-offset
     * `off` (2=t, 3=τ) crosses an integer multiple of dv. */
    private ticksFor(a: number[], off: number, dv: number, triangle: boolean): void {
        const ctx = this.ctx;
        const S = TRAIL_STRIDE;
        const n = a.length / S;
        ctx.lineWidth = 1;
        ctx.strokeStyle = '#15151a';                       // dark outline for contrast
        ctx.fillStyle = triangle ? '#ff5db1' : '#ffd24a';  // τ = magenta, t = amber
        let m = Math.floor(a[off] / dv) + 1;
        for (let k = 1; k < n; k++) {
            const v0 = a[(k - 1) * S + off], v1 = a[k * S + off];
            const dvSeg = v1 - v0;
            while (m * dv <= v1) {
                const target = m * dv;
                if (target >= v0) {
                    const f = dvSeg !== 0 ? (target - v0) / dvSeg : 0;
                    const sxv = a[(k - 1) * S] + f * (a[k * S] - a[(k - 1) * S]);
                    const syv = a[(k - 1) * S + 1] + f * (a[k * S + 1] - a[(k - 1) * S + 1]);
                    const px = this.sx(sxv), py = this.sy(syv);
                    ctx.beginPath();
                    if (triangle) {
                        ctx.moveTo(px, py - 4);
                        ctx.lineTo(px - 3.5, py + 3);
                        ctx.lineTo(px + 3.5, py + 3);
                        ctx.closePath();
                    } else {
                        ctx.arc(px, py, 3, 0, Math.PI * 2);
                    }
                    ctx.fill(); ctx.stroke();
                }
                m++;
            }
        }
    }

    /* A small two-hand clock: gray hand = coordinate time t, cyan hand = proper
     * time τ.  Both sweep with period CLOCK_PERIOD; the angular gap between them
     * grows as the accumulated dilation t-τ.  Digital readout: τ, Δ=t-τ, dτ/dt. */
    private drawClock(cx: number, cy: number, t: number, tau: number, rate: number): void {
        const ctx = this.ctx;
        const r = 13;
        const face = (ang: number, len: number, color: string, w: number) => {
            ctx.strokeStyle = color; ctx.lineWidth = w;
            ctx.beginPath();
            ctx.moveTo(cx, cy);
            ctx.lineTo(cx + len * Math.sin(ang), cy - len * Math.cos(ang));
            ctx.stroke();
        };
        ctx.beginPath(); ctx.arc(cx, cy, r, 0, Math.PI * 2);
        ctx.fillStyle = 'rgba(10,12,16,0.85)'; ctx.fill();
        ctx.strokeStyle = '#556'; ctx.lineWidth = 1.5; ctx.stroke();
        /* Lag wedge: the sector swept between the proper-time and coordinate-time
         * hands.  Its opening angle is the accumulated dilation t-τ -- a direct,
         * at-a-glance "how far has this clock fallen behind" gauge. */
        const angT = (2 * Math.PI * t) / CLOCK_PERIOD;
        const angTau = (2 * Math.PI * tau) / CLOCK_PERIOD;
        if (angT - angTau > 1e-3) {
            ctx.beginPath();
            ctx.moveTo(cx, cy);
            ctx.arc(cx, cy, r - 3, angTau - Math.PI / 2, angT - Math.PI / 2);
            ctx.closePath();
            ctx.fillStyle = 'rgba(255,93,177,0.22)'; ctx.fill();
        }
        face(angT, r - 3, '#aab', 1.5);      // coordinate
        face(angTau, r - 5, '#66ffcc', 2);   // proper
        ctx.fillStyle = '#66ffcc';
        ctx.beginPath(); ctx.arc(cx, cy, 1.6, 0, Math.PI * 2); ctx.fill();
        ctx.font = '9px ui-monospace, monospace'; ctx.textAlign = 'left';
        ctx.fillStyle = '#9ab';
        ctx.fillText(`τ${tau.toFixed(1)}`, cx + r + 3, cy - 2);
        ctx.fillStyle = '#c98';
        ctx.fillText(`Δ${(t - tau).toFixed(1)}`, cx + r + 3, cy + 8);
        if (Number.isFinite(rate)) {
            ctx.fillStyle = '#789';
            ctx.fillText(`${rate.toFixed(3)}`, cx + r + 3, cy + 18);
        }
    }

    /* Quiver: arrows on a `spacing`-cell grid showing the local field
     * direction + magnitude.  Lengths are normalized so the strongest sampled
     * arrow spans ~0.9 of the cell gap; alpha fades with relative magnitude so
     * weak regions recede.  Arrows are centered on their sample point.
     *
     * Performance: arrows are bucketed by relative magnitude into a few
     * opacity levels, each accumulated into one Path2D, so the whole quiver
     * costs ~2*BUCKETS canvas calls instead of one stroke+fill per arrow
     * (which is what made dense grids slow). */
    private drawVectors(v: VectorData, spacing: number): void {
        const ctx = this.ctx;
        const W = this.W, H = this.H;
        const step = Math.max(2, Math.round(spacing));
        const start = Math.floor(step / 2);

        /* First pass: peak magnitude among the sampled points. */
        let max = 0;
        for (let j = start; j < H; j += step) {
            for (let i = start; i < W; i += step) {
                const k = j * W + i;
                const m = Math.hypot(v.vx[k], v.vy[k]);
                if (m > max) max = m;
            }
        }
        if (max <= 0) return;

        const spacingPx = step * (this.cw / W);
        const scale = (spacingPx * 0.9) / max;   // px per field unit
        const head = Math.min(5, spacingPx * 0.22);

        const BUCKETS = 6;
        const lines: Path2D[] = [];
        const heads: Path2D[] = [];
        for (let b = 0; b < BUCKETS; b++) { lines.push(new Path2D()); heads.push(new Path2D()); }
        const invMax = 1 / max;

        for (let j = start; j < H; j += step) {
            for (let i = start; i < W; i += step) {
                const k = j * W + i;
                const fx = v.vx[k], fy = v.vy[k];
                const m = Math.hypot(fx, fy);
                if (m <= 0) continue;
                const rel = m * invMax;
                const b = Math.min(BUCKETS - 1, (rel * BUCKETS) | 0);
                const px = this.sx(i + 0.5), py = this.sy(j + 0.5);
                const dx = fx * scale, dy = -fy * scale;   // sim y is +up
                const x0 = px - dx * 0.5, y0 = py - dy * 0.5;
                const x1 = px + dx * 0.5, y1 = py + dy * 0.5;
                const lp = lines[b];
                lp.moveTo(x0, y0); lp.lineTo(x1, y1);
                const ang = Math.atan2(dy, dx);
                const hp = heads[b];
                hp.moveTo(x1, y1);
                hp.lineTo(x1 - head * Math.cos(ang - 0.4), y1 - head * Math.sin(ang - 0.4));
                hp.lineTo(x1 - head * Math.cos(ang + 0.4), y1 - head * Math.sin(ang + 0.4));
                hp.closePath();
            }
        }

        ctx.lineWidth = 1.25;
        for (let b = 0; b < BUCKETS; b++) {
            const a = (0.2 + 0.7 * ((b + 0.5) / BUCKETS)).toFixed(3);
            const color = `rgba(255,210,74,${a})`;
            ctx.strokeStyle = color; ctx.fillStyle = color;
            ctx.stroke(lines[b]); ctx.fill(heads[b]);
        }
    }

    private arrowHead(x: number, y: number, ang: number, color: string, s = 6): void {
        const ctx = this.ctx;
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.moveTo(x, y);
        ctx.lineTo(x - s * Math.cos(ang - 0.4), y - s * Math.sin(ang - 0.4));
        ctx.lineTo(x - s * Math.cos(ang + 0.4), y - s * Math.sin(ang + 0.4));
        ctx.closePath(); ctx.fill();
    }
}
