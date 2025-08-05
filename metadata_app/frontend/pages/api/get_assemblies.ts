import type { NextApiRequest, NextApiResponse } from "next";

export default async function handler(
  req: NextApiRequest,
  res: NextApiResponse
) {
  if (req.method !== "POST") {
    return res.status(405).json({ error: "Method not allowed" });
  }

  try {
    const payload = req.body;

    const response = await fetch("/api/assemblies/assemblies/filter", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      return res.status(500).json({ error: "Backend API failed" });
    }

    const data = await response.json();
    return res.status(200).json(data);
  } catch (error) {
    console.error("API error:", error);
    return res.status(500).json({ error: "Internal error" });
  }
}