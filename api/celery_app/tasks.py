from .celery_app import celery_app

@celery_app.task
def process_prompt(prompt: str):
    # Simulate some heavy computation
    return f"Processed prompt: {prompt}"
