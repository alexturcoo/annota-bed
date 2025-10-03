/** @type {import('next').NextConfig} */
const isDev = process.env.NODE_ENV !== 'production';

module.exports = {
  async rewrites() {
    if (isDev) {
      // Local dev: Next.js -> your local Flask (real annotator) on 5328
      return [{ source: '/api/:path*', destination: 'http://127.0.0.1:5328/:path*' }];
    }
    // Production: handled by Vercel serverless (demo mode)
    return [];
  },
};