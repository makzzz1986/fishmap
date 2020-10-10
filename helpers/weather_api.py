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

    def __init__(self, longitude, latitude, utc_timestamp):
        self.latitude = latitude
        self.longitude = longitude
        self.utc_timestamp = utc_timestamp

    def calculate_dang(self, height, period) -> int:
        self.dang = height * period * 5
        # print('Wave angle:', wave_spec['angle'], 'dang:', wave_spec['dang'])
        # print(wave_spec)
        return round(self.dang)

    # static setting to avoid consuming API calls
    def get_force(self, angle=45, height=0, period=0) -> dict:
        self.angle = angle
        self.height = height
        self.period = period
        return {'angle': round(self.angle), 'dang': self.calculate_dang(self.height, self.period)}

    def get_random(self, angle=45, height=0, period=0) -> dict:
        self.angle = angle + randint(-30, 30)
        self.height = randint(5, 30)/10
        self.period = randint(50, 100)/10
        return {'angle': round(self.angle), 'dang': self.calculate_dang(self.height, self.period)}

# Metod for [Stormglass API](https://stormglass.io/)
    def get_stormglass(self, api_token, utc) -> dict:
        response = requests.get(
            'https://api.stormglass.io/v2/weather/point',
            params={
                'lat': self.latitude,
                'lng': self.longitude,
                'params': ','.join([
                    'airTemperature',
                    'windDirection',
                    'windSpeed',
                    'waveDirection',
                    'waveHeight',
                    'wavePeriod',
                    'windWaveDirection',
                    'windWaveHeight',
                    'windWavePeriod'
                ]),
                'start': utc,  # Convert to UTC timestamp
                'end': utc  # we can request data for a few days, can we somehow use it to build few forecast at once?
            },
            headers={
                'Authorization': api_token
            }
        )

        json_data = response.json()
        print('Weather API response: ', json_data)
        self.angle = json_data['hours'][0]['waveDirection']['sg']
        return {'angle': round(self.angle), 'dang': self.calculate_dang(json_data['hours'][0]['waveHeight']['sg'], json_data['hours'][0]['wavePeriod']['sg'])}
