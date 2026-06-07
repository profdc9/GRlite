/* Lightweight hover tooltips for inspector field labels.
 *
 * Native `title` gives no control over the delay or styling, so this rolls a
 * tiny one: hover a label for DELAY ms and a small box appears beside it; the
 * next mouse movement (or leaving the label) dismisses it.  A single shared box
 * is reused for every label.  Small movements *while waiting* don't cancel the
 * timer (only leaving does), so you don't have to hold the cursor perfectly
 * still to summon it. */

const DELAY = 3000;   // ms to hover before the tip shows (per request: 3-4 s)

let box: HTMLDivElement | null = null;
let timer = 0;
let moveHandler: ((e: MouseEvent) => void) | null = null;

function getBox(): HTMLDivElement {
    if (!box) {
        box = document.createElement('div');
        box.className = 'tooltip-box';
        box.style.display = 'none';
        document.body.appendChild(box);
    }
    return box;
}

function clearTimer(): void { if (timer) { clearTimeout(timer); timer = 0; } }

function hide(): void {
    if (box) box.style.display = 'none';
    if (moveHandler) { document.removeEventListener('mousemove', moveHandler); moveHandler = null; }
}

function show(el: HTMLElement, text: string): void {
    const b = getBox();
    b.textContent = text;
    b.style.display = 'block';
    /* Position below the label, flipping above / clamping to stay on-screen. */
    const r = el.getBoundingClientRect();
    const bw = b.offsetWidth, bh = b.offsetHeight;
    let left = r.left;
    let top = r.bottom + 6;
    if (left + bw > window.innerWidth - 8) left = window.innerWidth - 8 - bw;
    if (top + bh > window.innerHeight - 8) top = r.top - 6 - bh;   // not enough room below
    b.style.left = `${Math.max(8, left)}px`;
    b.style.top = `${Math.max(8, top)}px`;
    /* Dismiss on the next real mouse movement. */
    moveHandler = () => hide();
    document.addEventListener('mousemove', moveHandler);
}

/* Attach a delayed hover tooltip to `el`.  No-op when text is empty. */
export function attachTooltip(el: HTMLElement, text: string | undefined): void {
    if (!text) return;
    el.style.cursor = 'help';
    el.addEventListener('mouseenter', () => {
        clearTimer();
        timer = window.setTimeout(() => { timer = 0; show(el, text); }, DELAY);
    });
    el.addEventListener('mouseleave', () => { clearTimer(); hide(); });
}
