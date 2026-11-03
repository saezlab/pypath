import { defineConfig } from 'eslint/config';
import js from '@eslint/js';
import ts from 'typescript-eslint';
import svelte from 'eslint-plugin-svelte';
import globals from 'globals';
import svelteConfig from './svelte.config.js';

export default defineConfig(
  {
    ignores: [
      '.svelte-kit/**',
      'build/**',
      'node_modules/**',
      '**/*.generated.ts',
      'playwright-report/**',
      'test-results/**',
    ],
  },
  js.configs.recommended,
  ts.configs.recommended,
  svelte.configs.recommended,
  {
    languageOptions: { globals: { ...globals.browser, ...globals.node } },
    rules: {
      '@typescript-eslint/no-unused-vars': [
        'error',
        { argsIgnorePattern: '^_', varsIgnorePattern: '^_', ignoreRestSiblings: true },
      ],
      '@typescript-eslint/no-explicit-any': 'off',
      'no-empty': ['error', { allowEmptyCatch: true }],
      // Maps/URLs are temporary calculations; this app deploys at the root path.
      // Keep these presentation preferences separate from correctness linting.
      'svelte/prefer-svelte-reactivity': 'off',
      'svelte/prefer-writable-derived': 'off',
      'svelte/no-navigation-without-resolve': 'off',
      'svelte/require-each-key': 'off',
    },
  },
  {
    files: ['src/lib/components/entity/MoleculeStructure.svelte'],
    rules: { 'svelte/no-dom-manipulating': 'off' },
  }, // OpenChemLib owns its SVG subtree.
  {
    files: ['**/*.svelte', '**/*.svelte.ts'],
    rules: { 'no-undef': 'off' },
    languageOptions: { parserOptions: { parser: ts.parser, svelteConfig } },
  },
);
