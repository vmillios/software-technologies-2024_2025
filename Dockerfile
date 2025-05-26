FROM python:3.12.10-slim

RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY ./requirements.txt .

RUN pip install -r requirements.txt

COPY . /app/

EXPOSE 8501

CMD [ "streamlit", "run", "./main.py", "--server.port=8501", "--server.address=0.0.0.0" ]