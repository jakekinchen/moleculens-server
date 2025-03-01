#!/usr/bin/env python3

"""
Free Throw Agentic System - Multi-Agent Version

A modular system for generating and visualizing animations with:
1. Static geometry code
2. Static animation code
3. LLM-generated captions
"""

import os
from typing import List, Dict, Any
from llm_service import LLMService, OpenAIProvider
from agents.geometry_agent import GeometryAgent
from agents.animation_agent import AnimationAgent
from agents.caption_agent import CaptionAgent
from agents.aggregator_agent import AggregatorAgent

def orchestrator(user_prompt: str) -> str:
    """
    1. Get geometry snippet (static)
    2. Get animation snippet (static)
    3. Get caption snippet (LLM-based)
    4. Merge all into a single HTML
    """
    print("\n[Orchestrator] Initializing LLM Service...")
    llm_service = LLMService(OpenAIProvider())

    geometry_agent = GeometryAgent(llm_service)
    animation_agent = AnimationAgent(llm_service)
    caption_agent = CaptionAgent(llm_service)
    aggregator_agent = AggregatorAgent()

    print("[Orchestrator] Retrieving geometry snippet...")
    geometry_snippet = geometry_agent.get_geometry_snippet(user_prompt)

    print("[Orchestrator] Retrieving animation snippet...")
    animation_snippet = animation_agent.get_animation_snippet(user_prompt)

    print("[Orchestrator] Generating caption snippet via LLM...")
    caption_snippet = caption_agent.get_caption_snippet(user_prompt)

    print("[Orchestrator] Combining into final HTML...")
    final_html = aggregator_agent.aggregate(
        geometry_snippet,
        animation_snippet,
        caption_snippet
    )

    # Write to output file
    output_file = "output.html"
    with open(output_file, "w") as f:
        f.write(final_html)
    print(f"\nVisualization saved to {output_file}")

    return final_html

def main():
    """Main entry point for the Free Throw Agentic System."""
    print("Welcome to the Multi-Agent Free Throw System!")
    print("This system combines static visualization code with LLM-generated captions.")
    
    # For demonstration, use a hardcoded prompt
    user_prompt = "animation of a magical reaction with multiple transitions"
    print(f"\nUsing prompt: {user_prompt}")

    try:
        final_code = orchestrator(user_prompt)
        print("\nOpen output.html in your browser to view the visualization.")
        
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
    except Exception as e:
        print(f"\nAn error occurred: {str(e)}")

if __name__ == "__main__":
    main()