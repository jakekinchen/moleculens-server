import json
import re
from datetime import timedelta

from bs4 import BeautifulSoup

# expects input like this: {"html": "<html>...</html>"}


def extract_captions_from_html(html_content):
    """Extract captions and their timecodes from HTML content.

    Returns a JSON structure with title and content array.
    """
    # Parse HTML content
    soup = BeautifulSoup(html_content, "html.parser")

    # Find all caption data in the JavaScript code
    script_content = soup.find_all("script")[-1].string  # Get the last script tag

    # Extract the captions array using regex
    captions_match = re.search(
        r"const\s+captions\s*=\s*\[(.*?)\];", script_content, re.DOTALL
    )
    if not captions_match:
        raise ValueError("Could not find captions array in the HTML")

    captions_text = captions_match.group(1)

    # Parse individual caption entries
    caption_entries = []
    for match in re.finditer(r"{[^}]*}", captions_text):
        entry_text = match.group(0)

        # Extract time and text using regex
        time_match = re.search(r"time:\s*(\d+)", entry_text)
        text_match = re.search(r'text:\s*"([^"]*)"', entry_text)

        if time_match and text_match:
            seconds = int(time_match.group(1))

            # Convert seconds to timecode format (HH:MM:SS)
            time_obj = timedelta(seconds=seconds)
            timecode = (
                f"{int(time_obj.total_seconds() // 3600):02d}:"
                f"{int((time_obj.total_seconds() % 3600) // 60):02d}:"
                f"{int(time_obj.total_seconds() % 60):02d}"
            )

            caption_entries.append(
                {
                    "timecode": timecode,
                    "description": f"Scene at {timecode}",
                    "caption": text_match.group(1),
                }
            )

    # Create the final JSON structure
    output = {
        "title": soup.title.string if soup.title else "Solar Wind Demonstration",
        "content": caption_entries,
    }

    return output


def convert_html_to_json(input_data, output_file_path=None):
    """Convert HTML content to JSON caption format.

    Input can be either a file path or a JSON object with an 'html' key.
    If output_file_path is provided, saves to file; otherwise returns
    JSON string.
    """
    try:
        # Check if input is a JSON object with HTML content
        if isinstance(input_data, dict) and "html" in input_data:
            html_content = input_data["html"]
        elif isinstance(input_data, str):
            # Check if input is a file path
            if input_data.endswith(".html"):
                with open(input_data, "r", encoding="utf-8") as file:
                    html_content = file.read()
            else:
                # Assume input is raw HTML string
                html_content = input_data
        else:
            raise ValueError(
                "Invalid input format. Expected file path, HTML string, or JSON object with 'html' key"
            )

        # Extract captions
        caption_data = extract_captions_from_html(html_content)

        # Convert to JSON string with proper formatting
        json_output = json.dumps(caption_data, indent=2, ensure_ascii=False)

        # Save to file if path provided
        if output_file_path:
            with open(output_file_path, "w", encoding="utf-8") as file:
                file.write(json_output)
            return f"JSON saved to {output_file_path}"

        return json_output

    except Exception as e:
        return f"Error processing HTML: {str(e)}"


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(
            "Usage: python html_to_caption_json.py <input_html_file_or_json> [output_json_file]"
        )
        sys.exit(1)

    # Check if input is JSON file
    try:
        with open(sys.argv[1], "r", encoding="utf-8") as f:
            input_data = json.load(f)
    except json.JSONDecodeError:
        # If not JSON, treat as HTML file
        input_data = sys.argv[1]
    except FileNotFoundError:
        print(f"Error: File {sys.argv[1]} not found")
        sys.exit(1)

    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    result = convert_html_to_json(input_data, output_file)
    print(result)
