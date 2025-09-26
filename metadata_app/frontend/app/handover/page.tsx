"use client";

import React, { useState, useEffect } from "react";
import { Accordion, AccordionItem, AccordionTrigger, AccordionContent } from "@/components/ui/accordion";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Terminal } from "lucide-react";
import { DataTable } from "@/app/tables/data-table-sorting";
import { Handover, columns } from "@/app/tables/handover_coumns";

export default function Page() {
  const [groupedData, setGroupedData] = useState<Record<string, Handover[]>>({});
  const [loading, setLoading] = useState(true);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);

  useEffect(() => {
  const fetchData = async () => {
    setLoading(true);
    setErrorMessage(null);

    try {
      const res = await fetch("/api/handover/handover/genebuilder");

      if (!res.ok) {
        const errorText = await res.text();
        throw new Error(`Server error: ${res.status} - ${errorText}`);
      }

      const json = await res.json();

      if (!json || typeof json !== "object") {
        throw new Error("API did not return a valid object");
      }

      setGroupedData(json); // API already grouped by genebuilder
    } catch (err) {
      setErrorMessage(err instanceof Error ? err.message : "Unknown error");
    } finally {
      setLoading(false);
    }
  };

  fetchData();
}, []);

  return (
    <div className="flex items-center justify-center">
      <div className="container mt-4 m-16 max-w-6xl">
        {errorMessage && (
          <Alert variant="destructive" className="mt-8">
            <Terminal />
            <AlertTitle>Error</AlertTitle>
            <AlertDescription>{errorMessage}</AlertDescription>
          </Alert>
        )}

        {loading ? (
          <p className="text-muted-foreground">Loading data...</p>
        ) : (
          <Accordion type="multiple" className="w-full">
            {Object.entries(groupedData).map(([person, data]) => (
              <AccordionItem key={person} value={person}>
                <AccordionTrigger>{person} ({data.length})</AccordionTrigger>
                <AccordionContent>
                  <div className="border-border border-2 rounded-md shadow-border mt-2">
                    <DataTable columns={columns} data={data} />
                  </div>
                </AccordionContent>
              </AccordionItem>
            ))}
          </Accordion>
        )}
      </div>
    </div>
  );
}