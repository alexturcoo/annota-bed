/** @type {import('next').NextConfig} */
const isDev = process.env.NODE_ENV !== 'production';

module.exports = {
  async rewrites() {
    if (isDev) {
      // Local dev: proxy Next.js -> local Flask (pnpm run flask-dev)
      return [{ source: '/api/:path*', destination: 'http://127.0.0.1:5328/:path*' }];
    }
    // Production on Vercel: /api/*.py runs as serverless. No rewrite.
    return [];
  },
};