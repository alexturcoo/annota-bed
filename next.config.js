/** @type {import('next').NextConfig} */
const isDev = process.env.NODE_ENV !== 'production';
const backendUrl = process.env.ANNOTATION_API_URL?.replace(/\/$/, '');

module.exports = {
  outputFileTracingRoot: __dirname,
  async rewrites() {
    if (backendUrl) {
      return [{ source: '/api/:path*', destination: `${backendUrl}/:path*` }];
    }
    if (isDev) {
      // Local dev: Next.js -> your local Flask (real annotator) on 5328
      return [{ source: '/api/:path*', destination: 'http://127.0.0.1:5328/:path*' }];
    }
    // Vercel's Python entry point receives the original /api/... URL.
    return [{ source: '/api/:path*', destination: '/api/index' }];
  },
};
