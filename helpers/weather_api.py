import requests
from random import randint


class Wave():
    latitude = 0.
    longitude = 0.
    angle = 0
    height = 0.
    period = 0.
    dang = 0.
    utc_timestamp = 0

    def __init__(self, latitude, longitude, utc_timestamp):
        self.longitude = longitude
        self.latitude = latitude
        self.utc_timestamp = utc_timestamp

    def calculate_dang(self, height, period) -> float:
        self.dang = height * period * 5
        # print('Wave angle:', wave_spec['angle'], 'dang:', wave_spec['dang'])
        # print(wave_spec)
        return self.dang

    # static setting to avoid consuming API calls
    def get_force(self, angle=45, height=0, period=0) -> float:
        self.angle = angle
        self.height = height
        self.period = period
        return self.calculate_dang(self.height, self.period)

    def get_random(self, angle=45, height=0, period=0) -> float:
        self.angle = angle + randint(-30, 30)
        self.height = randint(5, 30)/10
        self.period = randint(50, 100)/10
        return self.calculate_dang(self.height, self.period)

# Metod for [Stormglass API](https://stormglass.io/)

    def get_stormglass(self, api_token) -> float:
        response = requests.get(
        'https://api.stormglass.io/v2/weather/point',
        params={
            'lat': self.latitude,
            'lng': self.longitude,
            'params': ','.join([
                'airTemperature',
                'windDirection',
                'windSpeed',
                'swellDirection',
                'swellHeight',
                'swellPeriod',
                'secondarySwellPeriod',
                'secondarySwellDirection',
                'secondarySwellHeight',
                'waveDirection',
                'waveHeight',
                'wavePeriod',
                'windWaveDirection',
                'windWaveHeight',
                'windWavePeriod'
            ]),
            'start': self.utc_timestamp,  # Convert to UTC timestamp
            'end': self.utc_timestamp  # we can request data for a few days, can we somehow use it to build few forecast at once?
        },
        headers={
            'Authorization': api_token
        }
        )

        json_data = response.json()
        print(json_data)
