/* Minimal type declarations for the Emscripten-generated grlite.js module.
 * The .wasm is wrapped by an ES-module factory (MODULARIZE=1, EXPORT_ES6=1
 * in core/Makefile), so the default export is an async factory function. */

export interface GRliteModule {
  HEAPF32: Float32Array;
  HEAPU8: Uint8Array;
  HEAP32: Int32Array;
  _malloc(size: number): number;
  _free(ptr: number): void;
  cwrap<T extends (...args: unknown[]) => unknown>(
    name: string,
    returnType: string | null,
    argTypes: string[]
  ): T;
  ccall(
    name: string,
    returnType: string | null,
    argTypes: string[],
    args: unknown[]
  ): unknown;
  UTF8ToString(ptr: number): string;
  stringToUTF8(s: string, ptr: number, max: number): void;
}

declare const factory: () => Promise<GRliteModule>;
export default factory;
