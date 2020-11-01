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
        self.dang = height * period * 2
        # print('Wave angle:', wave_spec['angle'], 'dang:', wave_spec['dang'])
        # print(wave_spec)
        return round(self.dang)

    # static setting to avoid consuming API calls
    def get_force(self, angle=45, height=0, period=0) -> dict:
        self.angle = angle
        self.height = height
        self.period = period
        return {'angle': round(self.angle), 'dang': self.calculate_dang(self.height, self.period)}

    def get_random(self, angle=45, height=0, period=0, range_diff=20) -> dict:
        self.angle = angle + randint(range_diff*-1, range_diff)
        self.height = randint(5, 30)/10
        self.period = randint(50, 100)/10
        return {'angle': round(self.angle), 'dang': self.calculate_dang(self.height, self.period)}

# Metod for [Stormglass API](https://stormglass.io/)
    def get_stormglass(self, api_token, debug=False) -> dict:
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
                    'windWavePeriod',
                    'swellDirection',
                    'swellHeight',
                    'swellPeriod',
                    'secondarySwellPeriod',
                    'secondarySwellDirection',
                    'secondarySwellHeight'
                ]),
                'start': self.utc_timestamp,  # Convert to UTC timestamp
                'end': self.utc_timestamp  # we can request data for a few days, can we somehow use it to build few forecast at once?
            },
            headers={
                'Authorization': api_token
            }
        )

        json_data = response.json()

        if debug is True:
            print('Weather API response: ', json_data)
            
        api_angle = json_data['hours'][0]['waveDirection']['sg'] # 0 degree == 90 actual degree
        self.angle = (api_angle + 90) - ((api_angle + 90) // 360) * 360 # add 90 degree and remove 360 if it become more than 360
        return {'angle': round(self.angle), 
                'dang': self.calculate_dang(
                    json_data['hours'][0]['waveHeight']['sg'], 
                    json_data['hours'][0]['wavePeriod']['sg'])
                }
