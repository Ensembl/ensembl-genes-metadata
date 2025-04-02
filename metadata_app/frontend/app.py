import dash
from dash import html
import requests

app = dash.Dash(__name__)

# Fetch data from FastAPI backend
response = requests.get("http://backend:8000").json()

app.layout = html.Div([
    html.H1("Dash Frontend"),
    html.P(f"Backend says: {response['message']}")
])

if __name__ == "__main__":
    app.run_server(host="0.0.0.0", port=8050)