import type { NextConfig } from "next";
/** @type {import('next').NextConfig} */

const nextConfig: NextConfig = {
  eslint: {
    ignoreDuringBuilds: true, // Skip linting errors
  },
  typescript: {
    ignoreBuildErrors: true, // Skip type errors
  },
};

export default nextConfig;

