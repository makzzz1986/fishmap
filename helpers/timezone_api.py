import requests
from datetime import datetime, date, time, timedelta

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
    time_query = {}
    utc_offset = 0
    zone_name = ''
    country_name = ''

    def __init__(self, latitude, longitude, time_query):
        self.longitude = longitude
        self.latitude = latitude
        self.time_query = time_query # example: tomorrow at 8:00: {'days': 1, 'hours': 8}

    def get_static(self) -> dict:
        self.utc_offset = 3600
        future_timestamp = datetime.combine(date.today(), 
            datetime.min.time())+timedelta(days=self.time_query.get('days', 0), hours=self.time_query.get('hours', 0))
        return {'timestamp': int(future_timestamp.timestamp())+self.utc_offset, 'country_name': 'Portugal'}

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


a = Time(38.8, -9, {'days': 1, 'hours': 8})
print(a.get_static())