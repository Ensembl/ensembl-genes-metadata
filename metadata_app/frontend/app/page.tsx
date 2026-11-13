"use client";
import React from "react"
import {CardsStats} from "@/components/ui/cards_stats"
import {WelcomeCard} from "@/components/ui/card_welcome"
import {AnnotationGrowthChart} from "@/components/ui/card_annotation_growth"
import {AnnotationDurationChart} from "@/components/ui/card_annotation_duration"

export default function Page() {
    return (
        <div className="min-h-screen flex justify-center px-4 sm:px-6 lg:px-8">
            <div className="w-full max-w-6xl py-8 sm:py-16">
                <div className="pb-8">
                    <div className="w-full mb-8">
                        <WelcomeCard></WelcomeCard>
                    </div>
                    <div className="columns-1 space-y-4">
                        <div className="w-full mb-4 break-inside-avoid lg:break-before-column">
                            <CardsStats></CardsStats>
                        </div>
                        <div className="w-full mb-8">
                            <AnnotationDurationChart/>
                        </div>
                        <div className="w-full mb-8">
                            <AnnotationGrowthChart/>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    )
}
