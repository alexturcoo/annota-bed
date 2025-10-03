/** @type {import('next').NextConfig} */
const isDev = process.env.NODE_ENV !== 'production';

module.exports = {
  async rewrites() {
    if (isDev) {
      return [{ source: '/api/:path*', destination: 'http://127.0.0.1:5328/:path*' }];
    }
    return []; // prod: handled by Vercel functions
  },
};