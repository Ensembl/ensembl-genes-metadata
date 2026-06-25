"use client";
import React from "react"
import { CardsStats } from "@/components/ui/cards_stats"
import  { CardsDataTable } from "@/components/ui/card_projects"
import { WelcomeCard } from "@/components/ui/card_welcome"
import {ReportCard} from "@/components/ui/card_action_report";
import {DatabaseCleanup} from "@/components/ui/card_clean_up"
import {CardsDataTableGroup} from "@/components/ui/card_groups"

export default function Page() {
  return (
      <div className="min-h-screen flex justify-center px-4 sm:px-6 lg:px-8">
          <div className="w-full max-w-6xl py-10 sm:py-16">
            <div className="pt-4 pb-10 sm:pt-8">
                <div className="w-full mb-8">
                    <WelcomeCard></WelcomeCard>
                  </div>
              <div className="columns-1 space-y-4 lg:columns-2">

                  <div className="w-full mb-4 break-inside-avoid">
                    <CardsDataTable></CardsDataTable>
                  </div>
                  <div className="w-full mb-4 break-inside-avoid">
                    <DatabaseCleanup></DatabaseCleanup>
                  </div>
                  <div className="w-full mb-4 break-inside-avoid lg:break-before-column">
                    <CardsStats></CardsStats>
                  </div>
                 <div className="w-full mb-4 break-inside-avoid">
                    <CardsDataTableGroup></CardsDataTableGroup>
                 </div>
                  <div className="w-full mb-4 break-inside-avoid">
                    <ReportCard></ReportCard>
                  </div>
              </div>
            </div>
          </div>
      </div>
  )
}
