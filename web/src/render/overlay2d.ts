/* Canvas2D overlay drawn on top of the WebGL field canvas: particle trails,
 * markers, velocity arrows, and selection highlight.  Sim coordinates have
 * y=0 at the bottom (matching the WebGL particle pass); Canvas2D is y-down,
 * so we flip y in the mapping. */

import type { Particle } from '../sim/world';
import type { VectorData } from '../sim/fieldViews';

const VEL_PX_PER_C = 420;   // arrow pixels at |v| = c
const TRAIL_MAX = 600;       // points kept per particle

export interface Trails {
    /* trails[i] is a flat ring of [x0,y0,x1,y1,...] for particle i. */
    pts: number[][];
    push(particles: Particle[]): void;
    reset(n: number): void;
}

export function createTrails(): Trails {
    const t: Trails = {
        pts: [],
        reset(n: number) { t.pts = Array.from({ length: n }, () => []); },
        push(particles: Particle[]) {
            if (t.pts.length !== particles.length) t.reset(particles.length);
            for (let i = 0; i < particles.length; i++) {
                const arr = t.pts[i];
                arr.push(particles[i].x, particles[i].y);
                if (arr.length > TRAIL_MAX * 2) arr.splice(0, arr.length - TRAIL_MAX * 2);
            }
        },
    };
    return t;
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
           opts: { showTrails: boolean; showVelocity: boolean; selected: number;
                   vectors?: { data: VectorData; spacing: number } | null }): void {
        const ctx = this.ctx;
        ctx.clearRect(0, 0, this.cw, this.ch);

        /* Vector-field arrows (under the trails/particles). */
        if (opts.vectors) this.drawVectors(opts.vectors.data, opts.vectors.spacing);

        /* Trails: fade with age (oldest transparent). */
        if (opts.showTrails && trails) {
            ctx.lineWidth = 1.5;
            for (let i = 0; i < trails.pts.length; i++) {
                const a = trails.pts[i];
                const segs = a.length / 2;
                if (segs < 2) continue;
                for (let k = 1; k < segs; k++) {
                    const alpha = k / segs;
                    ctx.strokeStyle = `rgba(120,180,255,${(0.05 + 0.5 * alpha).toFixed(3)})`;
                    ctx.beginPath();
                    ctx.moveTo(this.sx(a[2 * (k - 1)]), this.sy(a[2 * (k - 1) + 1]));
                    ctx.lineTo(this.sx(a[2 * k]), this.sy(a[2 * k + 1]));
                    ctx.stroke();
                }
            }
        }

        for (const p of particles) {
            const x = this.sx(p.x), y = this.sy(p.y);

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
                ctx.strokeStyle = speed > 0.5 ? '#ff6060' : '#60a0ff';
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
            const sel = p.index === opts.selected;
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
        }
    }

    /* Quiver: arrows on a `spacing`-cell grid showing the local field
     * direction + magnitude.  Lengths are normalized so the strongest sampled
     * arrow spans ~0.9 of the cell gap; alpha fades with relative magnitude so
     * weak regions recede.  Arrows are centered on their sample point. */
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
        ctx.lineWidth = 1.25;

        for (let j = start; j < H; j += step) {
            for (let i = start; i < W; i += step) {
                const k = j * W + i;
                const fx = v.vx[k], fy = v.vy[k];
                const m = Math.hypot(fx, fy);
                if (m <= 0) continue;
                const px = this.sx(i + 0.5), py = this.sy(j + 0.5);
                const dx = fx * scale, dy = -fy * scale;   // sim y is +up
                const x0 = px - dx * 0.5, y0 = py - dy * 0.5;
                const x1 = px + dx * 0.5, y1 = py + dy * 0.5;
                const a = (0.2 + 0.7 * (m / max)).toFixed(3);
                ctx.strokeStyle = `rgba(255,210,74,${a})`;
                ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
                this.arrowHead(x1, y1, Math.atan2(y1 - y0, x1 - x0),
                               `rgba(255,210,74,${a})`, head);
            }
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
