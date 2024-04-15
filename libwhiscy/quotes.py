import random

_quotes = {"Mark Twain": "Too much of anything is bad, but too much good whiskey is barely enough.",
          "Raymond Chandler": "There is no bad whiskey. There are only some whiskeys that aren’t as good as others.",
          "Winston Churchill": "The water was not fit to drink. To make it palatable, we had to add whisky. By diligent effort, I learned to like it.",
          "Tommy Cooper": "I’m on a whisky diet. I’ve lost three days already.",
          "Humphrey Bogart": "I should never have switched from Scotch to Martinis.",
          "Johnny Carson": "Happiness is having a rare steak, a bottle of whisky, and a dog to eat the rare steak.",
          "Abraham Lincoln": "Tell me what brand of whiskey that Grant drinks. I would like to send a barrel of it to my other generals.",
          "Alexander Fleming": "A good gulp of hot whiskey at bedtime—it’s not very scientific, but it helps.",
          "W.C. Fields": "Always carry a flagon of whiskey in case of snakebite and furthermore always carry a small snake.",
          "James Joyce": "The light music of whiskey falling into a glass—an agreeable interlude.",
          "Igor Stravinsky": "My God, so much I like to drink Scotch that sometimes I think my name is Igor Stra-whiskey.",
          "Haruki Murakami": "The second whiskey is always my favorite. From the third on, it no longer has any taste. It's just something to pour into your stomach.",
          "Mikael Trellet": "Whiscy, to use without moderation.",
          "Jorge Roel": "Coffee keeps me going until it's acceptable to drink Whiskey"
          }


def get_one():
    one = random.choice(list(_quotes))
    return "\"{}\"\n  -  {}".format(_quotes[one], one)
