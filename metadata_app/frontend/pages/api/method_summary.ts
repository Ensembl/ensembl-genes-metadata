// pages/api/method_summary.ts
import type { NextApiRequest, NextApiResponse } from 'next'

export default async function handler(
  req: NextApiRequest,
  res: NextApiResponse
) {
  try {
    const response = await fetch("http://localhost:8000/api/home_page/home/method_summary")
    const data = await response.json()

    if (Array.isArray(data)) {
      res.status(200).json(data)
    } else {
      res.status(500).json({ error: "Data is not an array" })
    }
  } catch (error) {
    console.error("API proxy error:", error)
    res.status(500).json({ error: "Failed to fetch annotations data" })
  }
}