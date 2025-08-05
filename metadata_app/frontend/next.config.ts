import type { NextConfig } from "next";
/** @type {import('next').NextConfig} */

const nextConfig: NextConfig = {
  eslint: {
    ignoreDuringBuilds: true, // Skip linting errors
  },
  typescript: {
    ignoreBuildErrors: true, // Skip type errors
  },
  async rewrites() {
    return [
      {
        source: '/api/:path*',
        destination: 'http://127.0.0.1:8000/api/:path*', //  FastAPI server
      },
    ]
  },
};

export default nextConfig;

