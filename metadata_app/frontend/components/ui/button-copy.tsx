"use client";

import { useState } from "react";
import { motion, AnimatePresence } from "motion/react";
import { CheckIcon, CopyIcon } from "lucide-react";

import { Button } from "@/components/ui/button";
import { cn } from "@/lib/utils";

const MotionButton = motion.create(Button);
const MotionCopyIcon = motion.create(CopyIcon);
const MotionCheckIcon = motion.create(CheckIcon);

const copyIconVariants = {
  initial: { rotate: 0, scale: 1, y: 0 },
  hover: {
    rotate: [0, -8, 8, -4, 0],
    scale: 1.1,
    y: -0.5,
    transition: { duration: 0.4, ease: "easeInOut" },
  },
};

const checkIconVariants = {
  initial: { scale: 0, rotate: -30 },
  animate: {
    scale: [0, 1.3, 0.95, 1],
    rotate: 0,
    transition: { duration: 0.4, ease: "easeOut" },
  },
};

type CopyButtonProps = {
  text: string;
  className?: string;
  copiedLabel?: string;
  defaultLabel?: string;
  disabled?: boolean;
  onCopied?: () => void;
};

async function copyText(text: string) {
  try {
    if (navigator.clipboard?.writeText) {
      await navigator.clipboard.writeText(text);
      return true;
    }
  } catch (error) {
    console.error("Clipboard API write failed:", error);
  }

  const textarea = document.createElement("textarea");
  textarea.value = text;
  textarea.setAttribute("readonly", "");
  textarea.style.position = "fixed";
  textarea.style.left = "-9999px";
  textarea.style.top = "0";
  textarea.style.opacity = "0";
  document.body.appendChild(textarea);
  textarea.focus();
  textarea.select();
  textarea.setSelectionRange(0, textarea.value.length);

  let copied = false;
  try {
    copied = document.execCommand("copy");
  } finally {
    document.body.removeChild(textarea);
  }

  if (!copied) {
    window.prompt("Copy to clipboard:", text);
    return false;
  }

  return true;
}

export function CopyButton({
  text,
  className,
  copiedLabel = "Copied!",
  defaultLabel = "Copy",
  disabled = false,
  onCopied,
}: CopyButtonProps) {
  const [copied, setCopied] = useState(false);

  const handleCopy = async () => {
    if (!text.trim()) return;

    await copyText(text);
    setCopied(true);
    onCopied?.();
    window.setTimeout(() => setCopied(false), 1500);
  };

  return (
    <MotionButton
      type="button"
      variant="outline"
      className={cn(
        "relative cursor-pointer overflow-hidden disabled:opacity-100",
        className,
      )}
      onClick={handleCopy}
      disabled={disabled || copied}
      whileHover="hover"
      whileTap="tap"
      variants={{
        hover: { scale: 1.02 },
        tap: { scale: 0.98 },
      }}
      transition={{ type: "spring", stiffness: 400, damping: 25 }}
    >
      <span className="pointer-events-none invisible inline-flex items-center gap-1.5">
        <CopyIcon className="size-3.5" />
        {defaultLabel}
      </span>
      <AnimatePresence mode="wait" initial={false}>
        {copied ? (
          <motion.span
            key="copied"
            initial={{ opacity: 0, y: 6 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: -6 }}
            transition={{ duration: 0.15 }}
            className="absolute inset-0 inline-flex items-center justify-center gap-1.5 text-teal-400"
          >
            <MotionCheckIcon
              className="size-3.5 stroke-teal-400"
              variants={checkIconVariants}
              initial="initial"
              animate="animate"
            />
            {copiedLabel}
          </motion.span>
        ) : (
          <motion.span
            key="copy"
            initial={{ opacity: 0, y: -6 }}
            animate={{ opacity: 1, y: 0 }}
            exit={{ opacity: 0, y: 6 }}
            transition={{ duration: 0.15 }}
            className="absolute inset-0 inline-flex items-center justify-center gap-1.5"
          >
            <MotionCopyIcon
              className="size-3.5"
              variants={copyIconVariants}
              initial="initial"
            />
            {defaultLabel}
          </motion.span>
        )}
      </AnimatePresence>
    </MotionButton>
  );
}
