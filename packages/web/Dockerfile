# syntax=docker.io/docker/dockerfile:1

FROM node:22-alpine AS base

FROM base AS deps
RUN apk add --no-cache libc6-compat
WORKDIR /app/packages/web

RUN corepack enable

COPY packages/web/package.json packages/web/pnpm-lock.yaml* ./
RUN pnpm i --frozen-lockfile

FROM base AS builder
RUN corepack enable
WORKDIR /app
COPY --from=deps /app/packages/web/node_modules ./packages/web/node_modules
COPY packages/web ./packages/web

WORKDIR /app/packages/web

ARG PUBLIC_ENVIRONMENT
ARG PUBLIC_DOMAIN
ARG API_SERVICE_URL
ARG DOCKERIZED
ENV PUBLIC_ENVIRONMENT=$PUBLIC_ENVIRONMENT
ENV PUBLIC_DOMAIN=$PUBLIC_DOMAIN
ENV API_SERVICE_URL=$API_SERVICE_URL
ENV DOCKERIZED=$DOCKERIZED

RUN pnpm build
RUN pnpm prune --prod

FROM base AS runner
WORKDIR /app

ENV NODE_ENV=production

RUN addgroup --system --gid 1001 nodejs
RUN adduser --system --uid 1001 sveltekit

COPY --from=builder /app/packages/web/build ./build
COPY --from=builder /app/packages/web/node_modules ./node_modules
COPY --from=builder /app/packages/web/package.json ./

USER sveltekit

EXPOSE 8082

ENV BODY_SIZE_LIMIT=16M
ENV PORT=8082
ENV HOSTNAME="0.0.0.0"

HEALTHCHECK --interval=30s --timeout=10s --start-period=40s --retries=3 \
    CMD node -e "require('http').get('http://localhost:8082/api/health', (r) => r.statusCode === 200 ? process.exit(0) : process.exit(1))"

CMD ["node", "build"]