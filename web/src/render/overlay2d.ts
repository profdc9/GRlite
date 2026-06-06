/* Canvas2D overlay drawn on top of the WebGL field canvas: particle trails,
 * markers, velocity arrows, and selection highlight.  Sim coordinates have
 * y=0 at the bottom (matching the WebGL particle pass); Canvas2D is y-down,
 * so we flip y in the mapping. */

import type { Particle } from '../sim/world';

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
           opts: { showTrails: boolean; showVelocity: boolean; selected: number }): void {
        const ctx = this.ctx;
        ctx.clearRect(0, 0, this.cw, this.ch);

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

    private arrowHead(x: number, y: number, ang: number, color: string): void {
        const ctx = this.ctx, s = 6;
        ctx.fillStyle = color;
        ctx.beginPath();
        ctx.moveTo(x, y);
        ctx.lineTo(x - s * Math.cos(ang - 0.4), y - s * Math.sin(ang - 0.4));
        ctx.lineTo(x - s * Math.cos(ang + 0.4), y - s * Math.sin(ang + 0.4));
        ctx.closePath(); ctx.fill();
    }
}
