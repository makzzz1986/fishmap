import requests

# curl "https://api.timezonedb.com/v2.1/get-time-zone?key=&format=json&by=position&lat=38.0002&lng=-9"
#     {"status":"OK",
#     "message":"",
#     "countryCode":"PT",
#     "countryName":"Portugal",
#     "zoneName":"Europe\/Lisbon",
#     "abbreviation":"WEST",
#     "gmtOffset":3600,
#     "dst":"1",
#     "zoneStart":1585443600,
#     "zoneEnd":1603587600,
#     "nextAbbreviation":"WET",
#     "timestamp":1602258182,
#     "formatted":"2020-10-09 15:43:02"}


class Time():
    latitude = 0.
    longitude = 0.
    utc_offset = 0
    zone_name = ''
    country_name = ''

    def __init__(self, latitude, longitude):
        self.longitude = longitude
        self.latitude = latitude

    def get_static(self) -> dict:
        return {'utc_offset': 3600, 'country_name': 'Portugal'}

    def get_timezonedb(self, api_token) -> dict:
        response = requests.get(
            'https://api.timezonedb.com/v2.1/get-time-zone',
            params={
                'lat': self.latitude,
                'lng': self.longitude,
                'format': 'json',
                'key': api_token,
                'by': 'position'
            },
            headers={
                'Authorization': api_token
            }
        )

        json_data = response.json()
        print(json_data)
        self.utc_offset = json_data['gmtOffset']
        self.zone_name = json_data['zoneName']
        self.country_name = json_data['countryName']
        return {'utc_offset': self.utc_offset, 'country_name': self.country_name}
