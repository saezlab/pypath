import { API_SERVICE_URL } from '$lib/server/api-config';
import type { RequestHandler } from '@sveltejs/kit';
import { buildUpstreamUrl as upstreamUrl, proxyRequest } from './proxy-core';

export const buildUpstreamUrl = (pathname: string, search: string) =>
  upstreamUrl(pathname, search, API_SERVICE_URL);
export const forwardToBackend: RequestHandler = ({ request, url }) =>
  proxyRequest(request, url, API_SERVICE_URL);
export const forwardToAdminBackend = forwardToBackend;
