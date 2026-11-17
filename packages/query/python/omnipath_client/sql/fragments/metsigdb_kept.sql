-- A membership row {{input}} that the exclusion leaves in.
(op_unset({{exclude_classes}}) OR NOT list_has_any({{input}}.exclusion_classes, {{exclude_classes}}))
