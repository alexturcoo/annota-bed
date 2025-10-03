/** @type {import('next').NextConfig} */
const isDev = process.env.NODE_ENV !== "production";

module.exports = {
  async rewrites() {
    if (isDev) {
      // Local development: proxy to your locally running Flask
      return [{ source: "/api/:path*", destination: "http://127.0.0.1:5328/:path*" }];
    }
    // In production, let Vercel handle /api/* with serverless Python functions
    return [];
  },
};