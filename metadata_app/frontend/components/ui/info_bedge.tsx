import React from "react";
import { Badge } from "@/components/ui/badge";
import { CircleX } from "lucide-react";

interface MetricBadgeProps {
  label: string;
  onRemove: () => void;
}

export const MetricBadge: React.FC<MetricBadgeProps> = ({ label, onRemove }) => {
  return (
    <Badge
      variant="secondary"
      className="cursor-pointer group"
    >
      {label}
      <CircleX
        size={14}
        className="opacity-60 group-hover:opacity-100"
        onClick={(e) => {
          e.stopPropagation();
          onRemove();
        }}
      />
    </Badge>
  );
};