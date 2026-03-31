"use client";

import React, {useEffect, useState} from "react";

import {StatusBox, DataItem, PendingItem} from "@/components/ui/blocks/status_box"
import {HoReadyBox, Handover} from "@/components/ui/blocks/ho_ready_box"
import GridList02 from "@/components/ui/blocks/select_user"


export default function Page() {
    const title = "Genebuilder handover helper";
  const description =
    "This page lets you monitor handover ready cores and warns you about long-pending in progress annotations.";

  const [selectedGenebuilder, setSelectedGenebuilder] = useState<string | null>(null);

  const handleGenebuilderChange = (user: {
    name: string;
    role: string;
    imageUrl: string;
  }) => {
    setSelectedGenebuilder(user.role);
    // Save to browser
  localStorage.setItem("selectedGenebuilder", user.role);
  };
  const [hoTableData, setHoTable] = useState<Handover[]>([]);
    const [horeadyCount, setCountHOR] = useState<number>(0);
    const [dataCount, setCountData] = useState<number>(0);
    const [pendingCount, setPending] = useState<number>(0);
    const [listData, setListData] = useState<DataItem[]>([]);
    const [listPending, setListPending] = useState<PendingItem[]>([]);

  const [loading, setLoading] = useState(false);


  const handleGetHO = async () => {
    if (!selectedGenebuilder) return;
    setLoading(true);
    try {
      const payload = { genebuilder: selectedGenebuilder }; // string

      const res = await fetch("/api/handover/handover/genebuilder", {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          accept: "application/json",
        },
        body: JSON.stringify(payload),
      });

      const result = await res.json();
      console.log("API response:", result);

      if (result.df_ready) {
        setHoTable(result.df_ready);
        setCountHOR(result.count_ho_ready);
        setCountData(result.count_data);
        setPending(result.count_pending);
        setListData(result.list_data);
        setListPending(result.list_pending);


      } else {
        alert("No data found.");
      }

    } catch (error) {
      console.error("Error fetching data:", error);
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
  const saved = localStorage.getItem("selectedGenebuilder");
  if (saved) {
    setSelectedGenebuilder(saved);
  }
}, []);

  useEffect(() => {
    if (selectedGenebuilder) handleGetHO();
  }, [selectedGenebuilder]);


  return (
    <div className="flex items-center justify-center mt-15">
        <div>
        <div className="grid max-w-6xl gap-8 grid-cols-2">
        <h1 className="scroll-m-20 text-4xl font-extrabold tracking-tight text-balance">{title}</h1>
        <div className="max-w-lg min-w-lg flex justify-end">
            <GridList02 value={selectedGenebuilder} onSelect={handleGenebuilderChange} placeholder="Select a Genebuilder" />
            </div>
            </div>
            <p className="leading-7 [&:not(:first-child):mt-6]">{description}</p>
      {selectedGenebuilder ? (
            <div className="grid max-w-6xl gap-8">
        <StatusBox
          horeadyCount={horeadyCount ?? 0}
          dataCount={dataCount ?? 0}
          pendingCount={pendingCount ?? 0}
          listData={listData}
          listPending={listPending}
        />
        <HoReadyBox
          data={hoTableData}
          loading={loading}
          genebuilder={selectedGenebuilder}
        />
          </div>
              ) : (
              <div className="mt-8 w-full max-w-6xl h-64 flex items-center justify-center border-2 border-dashed text-gray-400">
          Select a Genebuilder to see the handover data
        </div>
      )}
    </div>
        </div>
  );
}