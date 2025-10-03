const isDev = process.env.NODE_ENV !== 'production';
const FLASK_URL = process.env.FLASK_URL; // set in Vercel env vars

module.exports = {
  async rewrites() {
    if (isDev) {
      return [{ source: '/api/:path*', destination: 'http://127.0.0.1:5328/:path*' }];
    }
    return [{ source: '/api/:path*', destination: `${FLASK_URL}/:path*` }];
  },
};