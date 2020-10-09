import requests

# curl "http://api.timezonedb.com/v2.1/get-time-zone?key=&format=json&by=position&lat=38.0002&lng=-9"

class Time():
    latitude = 0.
    longitude = 0.
    utc_offset = 0
    zoneStart = 0
    zoneEnd = 0
    zoneName = ''
    countryName = ''

    def __init__(self, latitude, longitude):
        self.longitude = longitude
        self.latitude = latitude

    def get_timezonedb(self, api_token, ) -> dict:
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
        print(json_data)
        self.angle = json_data['hours'][0]['waveDirection']['sg']
        return {'angle': round(self.angle), 'dang': self.calculate_dang(json_data['hours'][0]['waveHeight']['sg'], json_data['hours'][0]['wavePeriod']['sg'])}
