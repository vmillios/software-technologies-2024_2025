FROM python:3.12.10-slim

WORKDIR /app

COPY ./requirements.txt .

RUN pip install -r requirements.txt