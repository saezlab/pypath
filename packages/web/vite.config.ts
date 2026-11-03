import tailwindcss from '@tailwindcss/vite';
import { sveltekit } from '@sveltejs/kit/vite';
import { defineConfig } from 'vite';

// Use the same fixed-origin SvelteKit proxy in development and production.
export default defineConfig({ plugins: [tailwindcss(), sveltekit()] });
