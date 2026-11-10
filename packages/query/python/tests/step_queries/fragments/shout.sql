SELECT key, upper(label) || '!' AS shout FROM {{input}} WHERE _eid <= {{limit}}
