"""Extension module for LLM service to add image processing capabilities."""

from typing import Optional

from pydantic import BaseModel


class ImageToTextRequest(BaseModel):
    """Request structure for image-to-text processing."""

    image_data: str  # Base64 encoded image data
    prompt: Optional[str] = None
    max_tokens: Optional[int] = None


def extend_llm_service(module):
    """Extend the LLM service module with image processing capabilities."""
    # Add ImageToTextRequest to the module
    module.ImageToTextRequest = ImageToTextRequest

    # Get the original LLMService class
    original_llm_service = module.LLMService

    class ExtendedLLMService(original_llm_service):
        """Extended LLM service with image processing capabilities."""

        def process_image(self, request: ImageToTextRequest) -> str:
            """Process an image and return a text description."""
            # For now, we only support OpenAI's GPT-4 Vision model
            if self.config.provider != module.ProviderType.OPENAI:
                raise ValueError(
                    "Image processing is only supported with OpenAI provider"
                )

            # Use the centralized OpenAI client
            from ..utils.openai_client import get_client

            client = get_client(self.config.api_key)

            response = client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[
                    {
                        "role": "user",
                        "content": [
                            {
                                "type": "text",
                                "text": request.prompt or "What's in this image?",
                            },
                            {
                                "type": "image_url",
                                "image_url": {
                                    "url": f"data:image/jpeg;base64,{request.image_data}"
                                },
                            },
                        ],
                    }
                ],
                max_tokens=request.max_tokens or 300,
            )

            return response.choices[0].message.content or ""

    # Replace the original LLMService with the extended version
    module.LLMService = ExtendedLLMService
