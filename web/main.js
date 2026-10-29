// Loads the WebAssembly module and starts the simulator in #game.
import init, { start } from './pkg/nbody_gui.js';

const container = document.getElementById('game');

try {
  await init();
  await start(container);
} catch (err) {
  console.error(err);
  if (!container.querySelector('canvas')) {
    container.textContent = 'The simulator could not start in this browser: it needs WebAssembly and WebGL2.';
  }
}
